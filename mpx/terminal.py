import typer 
from pathlib import Path
from mpx.report import quality_control_consensus

app = typer.Typer()


@app.command()
def report(
    consensus_results: Path = typer.Argument(
        ..., help="Output path of the Nextflow pipeline for consensus genome assembly"
    )
):
    """
    Main report interface to generate the report file from one or multiple samples
    """
    quality_control_consensus(consensus_results=consensus_results)

app()