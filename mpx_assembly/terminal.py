import typer 
from pathlib import Path
from .report import quality_control_consensus

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
    Main report interface to generate the report file from one or multiple samples
    """
    quality_control_consensus(consensus_results=results)
