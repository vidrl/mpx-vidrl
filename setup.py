from setuptools import setup, find_packages

setup(
    name="mpx-assembly",
    url="https://github.com/vidrl/mpx-vidrl",
    author="Eike Steinig",
    author_email="eike.steinig@unimelb.edu.au",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "typer",
        "pyfastx",
        "rich",
        "pandas",
        "numpy",
        "seaborn",
        "matplotlib",
        "biopython",
        "cyvcf2"
    ],
    entry_points="""
    [console_scripts]
    mpx=report.terminal:app
    """,
    version="0.2.0",
    license="MIT",
    description="Monkeypox consensus assembly TWIST and ONT",
)