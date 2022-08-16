from setuptools import setup, find_packages

setup(
    name="mpx-assembly",
    url="https://github.com/esteinig/mpx-assembly",
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
        "matplotlib"
    ],
    entry_points="""
    [console_scripts]
    mpx=mpx_assembly.terminal:app
    """,
    version="0.1.0",
    license="MIT",
    description="Monkeypox assembly and all kinds of stuff",
)