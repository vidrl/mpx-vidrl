import typer 
from pathlib import Path
from .report import quality_control_consensus, variant_table

app = typer.Typer(add_completion=False)

report = typer.Typer()
app.add_typer(report, name="report")


@report.command()
def quality_control(
    results: Path = typer.Argument(
        ..., help="Output path of the Nextflow pipeline for consensus genome assembly"
    ),
    output: Path = typer.Option(
        "qc_table.tsv", help="QC table output file"
    ),
    subdir: str = typer.Option(
        None, help="Consensus sub directory"
    ),
    ont: bool = typer.Option(
        False, help="Output directory is from ARTIC"
    ),
    ont_ivar: bool = typer.Option(
        False, help="Output directory is from ONT IVAR workflow"
    ),
):
    """
    Quality control data
    """
    quality_control_consensus(results=results, subdir=subdir, table_output=output, ont=ont, ont_ivar=ont_ivar)


@report.command()
def variants(
        results: Path = typer.Argument(
            ..., help="Output path of the Nextflow pipeline for consensus genome assembly"
        ),
        output: Path = typer.Option(
            "qc_table.tsv", help="QC table output file"
        ),
        subdir: str = typer.Option(
            "high_freq", help="Consensus sub directory: low_freq | high_freq"
        ),
        genbank: Path = typer.Option(
            None, help="Genbank file of reference to use for additional variant annotations"
        ),
        min_completeness: float = typer.Option(
            95.0, help="Minimum completeness threshold to consider QC PASS"
        ),
        min_depth: float = typer.Option(
            0.0, help="Minimum depth threshold to consider QC PASS"
        ),
        mask: Path = typer.Option(
            None, help="Mask file (TSV) to annotate variants with masked regions"
        ),
        freq_alpha: float = typer.Option(
            0.3, help="Opacity of variant distribution dots"
        ),
        variant_pass: bool = typer.Option(
            False, help="Include only PASS variants in the summary outputs [mainly for low frequency summaries]"
        ),
        low_freq_depth: str = typer.Option(
            None, help="Dash separated string of {var_freq}:{min_depth}-{var_freq}:{min_depth} to filter variants below"
                       " frequency thresholds dependent on their depth"
        ),

):
    """
    Variant call summary
    """

    variant_table(
        results=results, subdir=subdir, variant_pass=variant_pass, low_freq_depth=low_freq_depth,
        min_complete=min_completeness/100, min_depth=min_depth, genbank_file=genbank, mask_file=mask,
        freq_alpha=freq_alpha
    )

app()