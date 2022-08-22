import typer 
from pathlib import Path
from .report import quality_control_consensus, snp_distance, variant_table

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
        "high_freq", help="Consensus sub directory: low_freq | high_freq"
    )
):
    """
    QC
    """
    quality_control_consensus(results=results, subdir=subdir, table_output=output)


@report.command()
def snp_dist(
    dist: Path = typer.Argument(
        ..., help="Pairwise SNP distance matrix in long form e.g. from PSDM"
    )
):
    """
    SNP distances
    """
    snp_distance(dist=dist)


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
        )
):
    """
    Variant call summary
    """

    variant_table(results=results, subdir=subdir, min_complete=0.95, min_depth=50, genbank_file=genbank)

