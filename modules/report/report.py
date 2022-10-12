""""
ARTIC-VIDRL report format

This report is mainly to summarize reference genome coverage and
identify primer issues from the coverage mask file for development.

Author: @esteinig
Date: 19/06/2022

"""
import dataclasses
from matplotlib.font_manager import json_load
import typer

from datetime import datetime
import json
import pandas
from numpy import histogram
from pathlib import Path
from typing import Optional, List, Tuple, Dict

from bokeh.plotting import figure
from bokeh.embed import components
from bokeh.palettes import Spectral6
from bokeh.models import ColumnDataSource, TableColumn, DataTable, HoverTool
from jinja2 import FileSystemLoader, Environment

@dataclasses.dataclass
class Primer:
    reference: str
    start: int
    end: int
    name: str
    pool: int
    strand: str
    sequence: Optional[str]


@dataclasses.dataclass
class NanoqReport:
    reads: int
    bases: int
    n50: int
    longest: int
    shortest: int
    mean_length: int
    median_length: int
    mean_quality: float
    median_quality: float
    length_thresholds: dict
    quality_thresholds: dict
    top_lengths: List[int]
    top_qualities: List[float]
    filtered: int


@dataclasses.dataclass
class Mask:
    reference: str
    start: int
    end: int


@dataclasses.dataclass
class MaskedPrimers:
    mask: Mask
    primers: List[Primer]


@dataclasses.dataclass
class SampleData:
    sample_id: str
    coverage_file: Path
    mask_file: Optional[Path]
    nanoq_file: Path
    length_file: Path
    quality_file: Path


def collect_samples(coverage_files: List[Path], mask_files: Optional[List[Path]], nanoq_files: List[Path], read_length_files: List[Path], read_qualities_files: List[Path]) -> List[SampleData]:
    """
    Collect samples by sample name from matching mask and coverage files
    """

    # Sample names should conform to pipeline outputs being
    # the prefix in a '.' delimited name, we ae using the
    # coverage file name as it should always be present

    samples = []
    for cov_file in coverage_files:
        try:
            sample_id = cov_file.name.split('.')[0]
        except IndexError:
            raise ValueError(f"Could not find sample name in file: {cov_file}")

        if mask_files:
            try:
                mask_file = [
                    mask_file for mask_file in mask_files
                    if mask_file.name.startswith(sample_id)
                ][0]
            except IndexError:
                print(f"Could not find matching mask file for sample: {sample_id}")
                mask_file = None
        else:
            mask_file = None


        try:
            nanoq_file = [
                nanoq_file for nanoq_file in nanoq_files
                if nanoq_file.name.startswith(sample_id)
            ][0]
        except IndexError:
            print(f"Could not find matching nanoq file for sample: {sample_id}")
            nanoq_file = None

        try:
            length_file = [
                length_file for length_file in read_length_files
                if length_file.name.startswith(sample_id)
            ][0]
        except IndexError:
            print(f"Could not find matching read length file for sample: {sample_id}")
            length_file = None

        try:
            quality_file = [
                quality_file for quality_file in read_qualities_files
                if quality_file.name.startswith(sample_id)
            ][0]
        except IndexError:
            print(f"Could not find matching read quality file for sample: {sample_id}")
            quality_file = None

        samples.append(
            SampleData(sample_id=sample_id, coverage_file=cov_file, mask_file=mask_file, nanoq_file=nanoq_file, length_file=length_file, quality_file=quality_file)
        )
    return samples


def read_scheme_bed(file: Path, sep: str = "\t") -> pandas.DataFrame:

    scheme = pandas.read_csv(file, sep=sep, header=None)

    scheme_column_names = ["reference", "start", "end", "name", "pool", "strand"]
    if len(scheme.columns) == 6:
        scheme.columns = scheme_column_names
    elif len(scheme.columns) == 7:
        scheme.columns = scheme_column_names + ["sequence"]  # internal development format
    else:
        raise ValueError("Scheme file requires six or seven columns.")

    return scheme


def read_coverage_mask(file: Path) -> pandas.DataFrame:
    """
    Assuming 0-based, half-open format as in BED (does not include end position)
    """
    return pandas.read_csv(
        file, sep="\t", header=None, names=["reference", "start", "end"],
    )


def read_coverage_bed(file: Path, sep: str = "\t") -> pandas.DataFrame:
    """
    Read output from `covtobed`
    """
    return pandas.read_csv(file, sep=sep, header=None, names=["reference", "start", "end", "coverage"])


def get_affected_primers(mask_txt: Path, scheme_bed: Path) -> List[MaskedPrimers]:
    """
    Compare the `coverage_mask.txt` output file with the `.scheme.bed` file
    and report specific primers involved in low coverage masked sites in the
    final assembly.
    """

    primer_scheme = read_scheme_bed(file=scheme_bed)
    coverage_mask = read_coverage_mask(file=mask_txt)

    affected = []
    for _, mask in coverage_mask.iterrows():
        mask_reference = mask["reference"]
        mask_start, mask_end = mask['start'], mask['end']
        mask_range = range(mask_start, mask_end)

        # Due to multiple chromosome format of BED file subset
        # the scheme data to the correct reference sequence name
        primer_scheme_ref = primer_scheme[primer_scheme["reference"] == mask_reference]
        if primer_scheme_ref.empty:
            raise ValueError(f"Could not find {mask_reference} in primer scheme file")

        _mask = Mask(
            reference=mask_reference,
            start=mask_start,
            end=mask_end
        )

        affected_primers = []
        _seen = []
        for _, primer in primer_scheme_ref.iterrows():
            primer_range = range(primer['start'], primer['end'])
            primer_dict = primer.to_dict()
            primer_name = primer_dict["name"]

            for position in mask_range:
                if position in primer_range:
                    if primer_name not in _seen:
                        affected_primers.append(
                            Primer(
                                reference=primer_dict["reference"],
                                start=primer_dict["start"],
                                end=primer_dict["end"],
                                name=primer_dict["name"],
                                pool=primer_dict["pool"],
                                strand=primer_dict["strand"],
                                sequence=primer.get("sequence", None)
                            )
                        )
                        _seen.append(primer_name)

        affected.append(MaskedPrimers(mask=_mask, primers=affected_primers))

    return affected


def get_masked_primer_data(masked_primers: List[MaskedPrimers]) -> pandas.DataFrame:
    """
    Get primers involved in mask as rows for a table in Bokeh
    """
    data = []
    for masked_primer in masked_primers:
        for primer in masked_primer.primers:
            data.append(
                {
                    'reference': masked_primer.mask.reference,
                    'masked_region': f"{masked_primer.mask.start} - {masked_primer.mask.end}",
                    'masked_bp': f"{masked_primer.mask.end - masked_primer.mask.start}",
                    'primer': primer.name,
                    'region': f"{primer.start} - {primer.end}",
                    'pool': primer.pool
                }
            )

    return pandas.DataFrame(data)


def get_primer_pools(primer_scheme: pandas.DataFrame) -> Dict[int, List[Primer]]:
    """
    Get the primers in a dict by pool [key]:
    """
    pool_primers = {}
    for _, primer in primer_scheme.iterrows():
        primer_dict = primer.to_dict()
        primer = Primer(
            reference=primer_dict["reference"],
            start=primer_dict["start"],
            end=primer_dict["end"],
            name=primer_dict["name"],
            pool=primer_dict["pool"],
            strand=primer_dict["strand"],
            sequence=primer.get("sequence", None)
        )
        if primer.pool not in pool_primers.keys():
            pool_primers[primer.pool] = [primer]
        else:
            pool_primers[primer.pool].append(primer)

    return pool_primers


def get_sample_coverage_data(cov_bed: Path, scheme_bed: Path) -> Tuple[List, List, List]:
    """
    Get coverage data from `covtobed` output and transform for plotting with Bokeh
    """
    cov = read_coverage_bed(file=cov_bed)
    primer_scheme = read_scheme_bed(file=scheme_bed)

    # Get the primer scheme ranges by pool
    pool_primers = get_primer_pools(primer_scheme=primer_scheme)

    full_x = []
    full_y = []
    full_pool = []
    for _, region in cov.iterrows():
        for i in range(region['start'], region['end']):
            full_x.append(i)
            # Check which primer pool the base is in and add
            # it to the list of pool values (this is really
            # inefficient right now)
            base_pool = None
            for pool, primers in pool_primers.items():
                for primer in primers:
                    if i in range(primer.start, primer.end):
                        base_pool = pool

            if base_pool:
                base_pool = str(base_pool)
            full_pool.append(base_pool)
            full_y.append(int(region['coverage']))

    assert len(full_x) == len(full_pool)

    return full_x, full_y, full_pool


def get_numeric(sample_id: str):
    """
    If numerics in sample id get those to sort, otherwise return string
    """
    digits = [c for c in sample_id if c.isdigit()]
    if digits:
        try:
            return int("".join(digits))
        except ValueError:
            raise ValueError(f"Could not convert digits in sample name: {digits}")
    else:
        return 0


def read_params(params_file: Path) -> dict:

    with params_file.open() as pf:
        pipeline_params = json.load(pf)

    return pipeline_params


def read_txt_file(file: Path) -> List:

    with file.open() as infile:
        data = [line.strip() for line in infile]
    return data


def read_nanoq_json(file: Path) -> NanoqReport:
    """
    Read the nanoq report file into dataclass
    """
    with file.open() as infile:
        nanoq_dict = json.load(infile)
    return NanoqReport(
        reads=nanoq_dict['reads'],
        bases=nanoq_dict['bases'],
        n50=nanoq_dict['n50'],
        longest=nanoq_dict['longest'],
        shortest=nanoq_dict['shortest'],
        mean_length=nanoq_dict['mean_length'],
        median_length=nanoq_dict['median_length'],
        mean_quality=nanoq_dict['mean_quality'],
        median_quality=nanoq_dict['median_quality'],
        length_thresholds=nanoq_dict['length_thresholds'],
        quality_thresholds=nanoq_dict['quality_thresholds'],
        top_lengths=nanoq_dict['top_lengths'],
        top_qualities=nanoq_dict['top_qualities'],
        filtered=nanoq_dict['filtered']
    )


def get_read_length_plot(read_lengths: List[int]):


    read_lengths = [int(rl) for rl in read_lengths]
    arr_hist, edges = histogram(read_lengths, bins=20)

    # Column data source
    arr_df = pandas.DataFrame({'count': arr_hist, 'left': edges[:-1], 'right': edges[1:]})
    arr_df['f_count'] = ['%d' % count for count in arr_df['count']]
    arr_df['f_interval'] = ['%d to %d ' % (left, right) for left, right in zip(arr_df['left'], arr_df['right'])]

    # column data source
    arr_src = ColumnDataSource(arr_df)

    # Set up the figure same as before
    p = figure(plot_width=320,
               plot_height=320,
               tools="pan,wheel_zoom,save,reset",
               title = "Read length distribution (all)",
               x_axis_label = "bp")

    # Add a quad glyph with source this time
    p.quad(
        bottom=0, 
        top='count', 
        left='left', 
        right='right', 
        source=arr_src,
        fill_color=Spectral6[0],
        hover_fill_alpha=0.7,
        hover_fill_color=Spectral6[0],
        line_color='black'
    )

    # Add style to the plot
    p.title.align = 'center'

    # Add a hover tool referring to the formatted columns
    hover = HoverTool(tooltips = [("Read length", '@f_interval'),
                                  ('Count', '@f_count')])

    # Add the hover tool to the graph
    p.add_tools(hover)

    return p


def get_read_qualities_plot(read_qualities: List[float]):

    read_quals = [float(rl) for rl in read_qualities]
    arr_hist, edges = histogram(read_quals, bins=20)

    # Column data source
    arr_df = pandas.DataFrame({'count': arr_hist, 'left': edges[:-1], 'right': edges[1:]})
    arr_df['f_count'] = ['%d' % count for count in arr_df['count']]
    arr_df['f_interval'] = ['%d to %d ' % (left, right) for left, right in zip(arr_df['left'], arr_df['right'])]

    # column data source
    arr_src = ColumnDataSource(arr_df)

    # Set up the figure same as before
    p = figure(plot_width = 320, 
               plot_height = 320,
               tools="pan,wheel_zoom,save,reset",
               title = "Read quality distribution (all)",
               x_axis_label = "Q")

    # Add a quad glyph with source this time
    p.quad(
        bottom=0, 
        top='count', 
        left='left', 
        right='right', 
        source=arr_src,
        fill_color=Spectral6[0],
        hover_fill_alpha=0.7,
        hover_fill_color=Spectral6[0],
        line_color='black'
    )

    # Add style to the plot
    p.title.align = 'center'

    # Add a hover tool referring to the formatted columns
    hover = HoverTool(tooltips = [("Read quality", '@f_interval'),
                                  ('Count', '@f_count')])

    # Add the hover tool to the graph
    p.add_tools(hover)

    return p


def generate_test_html(
    scheme: Path, params: Path,
    coverage_files: List[Path],
    mask_files: Optional[List[Path]],
    nanoq_report_files: List[Path],
    read_length_files: List[Path],
    read_qualities_files: List[Path],
    barcodes: bool,
    report_cov: int
):

    loader = FileSystemLoader(searchpath=Path(__file__).parent)
    env = Environment(loader=loader)

    template_file = "template.html"
    template = env.get_template(template_file)

    sample_data = collect_samples(
        coverage_files=coverage_files, 
        mask_files=mask_files, 
        nanoq_files=nanoq_report_files, 
        read_length_files=read_length_files,
        read_qualities_files=read_qualities_files
    )

    params_dict = read_params(params_file=params)

    artic_ended = datetime.now().strftime("%Y-%m-%d %H:%M")
    params_dict["finished"] = artic_ended
    params_dict["report_coverage"] = report_cov
    
    print(params_dict)

    # Get the primer scheme ranges by pool
    # pool_primers = get_primer_pools(primer_scheme=primer_scheme)

    samples = []

    if barcodes:
        sample_data = sorted(sample_data, key=lambda x: get_numeric(x.sample_id))

    total_read_lengths = []
    total_read_qualities = []
    for sample in sample_data:
        
        # Coverage plot
    
        x, y2, pool = get_sample_coverage_data(cov_bed=sample.coverage_file, scheme_bed=scheme)
        df = pandas.DataFrame(
            {'x': x, 'y': y2, 'pool': pool}
        )
        over_threshold = len([cov for cov in y2 if cov >= report_cov])/len(y2)
        p = figure(
            title=f"{sample.sample_id}: {over_threshold*100:.2f}% > {report_cov}x",
            x_axis_label="Position",
            y_axis_label="Depth",
            tools="pan,wheel_zoom,save,reset",
            width=500,
            height=450
        )
        source = ColumnDataSource(df)
        p.varea(x='x', y1=0, y2='y', source=source, fill_color=Spectral6[0], fill_alpha=1)
        plot_script, plot_div = components(p)

        # Affected primer table

        masked_primers = get_affected_primers(mask_txt=sample.mask_file, scheme_bed=scheme)
        masked_primer_data = get_masked_primer_data(masked_primers=masked_primers)
        
        columns = [
            TableColumn(field='masked_region', title='Masked region'),
            TableColumn(field='masked_bp', title='Bases'),
            TableColumn(field='primer', title='Affected primer'),
            TableColumn(field='pool', title='Pool'),
        ]
        table = DataTable(
            source=ColumnDataSource(masked_primer_data),
            columns=columns,
            autosize_mode='fit_columns'
        )

        table_script, table_div = components(table)

        # Read nanoq reports

        nanoq_sample_report = read_nanoq_json(file=sample.nanoq_file)

        # Read summary table
        stats_data = pandas.DataFrame(
            [
                ('Reads', str(nanoq_sample_report.reads)),
                ('Filtered', str(nanoq_sample_report.filtered)),
                ('Bases', str(nanoq_sample_report.bases)),
                ('Shortest', str(nanoq_sample_report.shortest)),
                ('Longest', str(nanoq_sample_report.longest)),
                ('Mean length', str(nanoq_sample_report.mean_length)),
                ('Median length', str(nanoq_sample_report.median_length)),
                ('Median quality', str(nanoq_sample_report.median_quality)),
                ('Mean quality', str(nanoq_sample_report.mean_quality)),
            ], columns=["stats", "value"]
        )

        print(stats_data)

        columns = [
            TableColumn(field='stats', title='Summary'),
            TableColumn(field='value', title='Value'),
        ]
        stats_table = DataTable(
            source=ColumnDataSource(stats_data),
            columns=columns,
            autosize_mode='fit_columns',
            index_position=None,
            width=200,
            height=250
        )

        stats_table_script, stats_table_div = components(stats_table)

        # Summarize to sample object for Jinja2 template

        samples.append(
            {
                'id': sample.sample_id,
                'plot_div': plot_div,
                'plot_script': plot_script,
                'table_div': table_div,
                'table_script': table_script,
                'stats_table_div': stats_table_div,
                'stats_table_script': stats_table_script,
            }
        )

        for read_length in read_txt_file(sample.length_file):
            total_read_lengths.append(read_length)

        for read_quality in read_txt_file(sample.quality_file):
            total_read_qualities.append(read_quality)
    
    p = get_read_length_plot(total_read_lengths)
    plot_length_script, plot_length_div = components(p)


    p = get_read_qualities_plot(total_read_qualities)
    plot_qual_script, plot_qual_div = components(p)

    # p = get_read_qualities_plot(total_read_qualities)
    # plot_quality_script, plot_quality_div = components(p)

    tmp = template.render(
        samples=samples,
        version=params_dict['version'],
        commit="7002f91",
        start_time=params_dict['started'],
        end_time=params_dict['finished'],
        n_samples=str(len(sample_data)),
        barcodes=params_dict['barcodes'],
        min_length=params_dict['min_length'],
        max_length=params_dict['max_length'],
        min_quality=params_dict['min_quality'],
        normalise=params_dict['normalise'],
        medaka_model=params_dict['medaka_model'],
        report_title=f": {params_dict['report_title']}",
        plot_length_div=plot_length_div,
        plot_length_script=plot_length_script,
        plot_qual_script=plot_qual_script,
        plot_qual_div=plot_qual_div,
        report_coverage=params_dict["report_coverage"]
    )

    with open("report.html", "w") as fh:
        fh.write(tmp)


app = typer.Typer()


@app.command()
def report(
    scheme: Path = typer.Option(..., help="Primer scheme BED file"),
    params: Path = typer.Option(None, help="Nextflow params file output from project."),
    glob_path: Path = typer.Option(Path.cwd(), help="The directory containing the report data files."),
    coverage_glob: str = typer.Option("*.coverage.bed", help="Coverage reference BED files"),
    mask_glob: str = typer.Option("*.coverage_mask.txt", help="Coverage mask TXT files"),
    read_lengths_glob: str = typer.Option("*.read_lengths.txt", help="Read lengths files from nanoq"),
    read_qualities_glob: str = typer.Option("*.read_qualities.txt", help="Read qualities files from nanoq"),
    nanoq_report_glob: str = typer.Option("*.nanoq.json", help="Nanoq sample report JSON files"),
    barcodes: bool = typer.Option(False, help="Barcode mode, sorts by trailing digits"),
    report_cov: int = typer.Option(20, help="Report coverage over this threshold per sample, average across reference"),
):
    """
    Main report interface to generate the report file from one or multiple samples
    """

    # Typer/Click do not support multiple inputs :(
    mask_files = [p for p in glob_path.glob(mask_glob)]
    coverage_files = [p for p in glob_path.glob(coverage_glob)]
    nanoq_report_files = [p for p in glob_path.glob(nanoq_report_glob)]
    read_lengths = [p for p in glob_path.glob(read_lengths_glob)]
    read_qualities = [p for p in glob_path.glob(read_qualities_glob)]


    generate_test_html(
        scheme=scheme,
        params=params, 
        coverage_files=coverage_files, 
        mask_files=mask_files,
        nanoq_report_files=nanoq_report_files,
        read_length_files=read_lengths,
        read_qualities_files=read_qualities,
        barcodes=barcodes,
        report_cov=report_cov
    )


if __name__ == "__main__":
    app()
