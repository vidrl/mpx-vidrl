import typer 
from pathlib import Path
from .report import quality_control_consensus, snp_distance

app = typer.Typer(add_completion=False)

report = typer.Typer()
app.add_typer(report, name="report")


@report.command()
def consensus(
    results: Path = typer.Argument(
        ..., help="Output path of the Nextflow pipeline for consensus genome assembly"
    )
):
    """
    QC
    """
    quality_control_consensus(consensus_results=results)


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
