import typer 
from pathlib import Path

app = typer.Typer()

@app.command()
def report(
    consensus: Path = typer.Argument(..., help="Output path of the Nextflow pipeline for consensus genome assembly")
):
    """
    Main report interface to generate the report file from one or multiple samples
    """
    
    pass
