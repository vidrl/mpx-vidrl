"""
Monkeypox assembly report
"""

from cmath import nan
import json
import numpy

import pandas
from pathlib import Path
from pyfastx import Fasta
from rich.table import Table
from rich import print as rprint
from dataclasses import dataclass
from typing import Optional, List, Tuple
import numpy as np
from statistics import median, mean
import seaborn as sns
from matplotlib import pyplot as plt


@dataclass
class SampleFiles:
    assembly: Path
    fastp: Path
    samtools: Path


@dataclass
class SampleQC:
    name: str
    reads: Optional[int]
    qc_reads: Optional[int]
    aligned_reads: Optional[int]
    coverage: Optional[float]
    mean_depth: Optional[float]
    missing_sites: Optional[int]
    completeness: Optional[float]

    def to_list(self):
        return [
            self.name,
            self.reads,
            self.qc_reads,
            self.aligned_reads,
            self.coverage,
            self.mean_depth,
            self.missing_sites,
            self.completeness
        ]


def get_fastp_data(file: Path) -> Tuple[int, int]:
    """
    Get fastp data - divive by two for paired-end reads
    """
    with file.open() as infile:
        fastp_data = json.load(infile)

    all_reads = fastp_data["summary"]["before_filtering"]["total_reads"]
    qc_reads = fastp_data["summary"]["after_filtering"]["total_reads"]

    return all_reads, qc_reads  # Illumina PE


def get_samtools_data(file: Path) -> Tuple[int, float, float]:
    """
    Get samtools coverage data
    """
    content = file.open().readlines()[1].strip().split("\t")
    return int(content[3]), round(float(content[5]), 4), round(float(content[6]), 6)  # numreads, coverage, meandepth


def get_consensus_assembly_data(file: Path) -> Tuple[float or None, int]:
    """
    Get consensus sequence and missing site proportion (N) - should only have a single sequence
    """

    seq_data = [seq for seq in Fasta(str(file), uppercase=True, build_index=False)]
    seq = seq_data[0][1]
    ncount = seq.count("N")
    try:
        completeness = round(100 - ((ncount / len(seq))*100), 6)
    except ZeroDivisionError:
        completeness = None

    return completeness, ncount


def create_rich_table(samples: List[SampleQC], title: str, patient_id: bool = True, table_output: Path = None):

    df = pandas.DataFrame(
        [sample.to_list() for sample in samples],
        columns=["Sample", "Reads", "QC Reads", "Alignments", "Coverage", "Mean Depth", "Missing", "Completeness"]
    )

    if not patient_id:
        df = df.sort_values(["Sample", "Completeness", "Coverage", "Mean Depth"])
    else:
        # Sort first by sample patient identifier then number of that patient sample
        # must comply with Mona's format: ID_{Patient}_{Number} e.g. MPX_A_1 and MPX_A_2
        patient_samples = {}
        for i, row in df.iterrows():
            sample_id = row["Sample"]
            sample_content = sample_id.split("_")
            patient_id = sample_content[1]
            sample_number = sample_content[2]

            if patient_id not in patient_samples.keys():
                patient_samples[patient_id] = [(int(sample_number), row.tolist())]
            else:
                patient_samples[patient_id].append((int(sample_number), row.tolist()))

        sorted_patient_samples = {}
        for patient_id, patient_data in patient_samples.items():
            sorted_patient_samples[patient_id] = sorted(patient_data, key=lambda x: x[0])

        sorted_samples = dict(sorted(sorted_patient_samples.items()))  # Python 3.7+

        df = pandas.DataFrame(
            [sample[1] for _, data in sorted_samples.items() for sample in data],
            columns=[
                "Sample",
                "Reads",
                "QC Reads",
                "Alignments",
                "Coverage",
                "Mean Depth",
                "Missing (N)",
                "Completeness"
            ]
        )

    table = Table(title=title)
    for cname in df.columns:
        if cname != "Sample":
            justify = "right"
        else:
            justify = "left"
        table.add_column(cname, justify=justify, no_wrap=False)
    for _, row in df.iterrows():
        if row["Completeness"] >= 99.9:
            row_color = "#66bd63"
        elif 95.0 <= row["Completeness"] < 99.9:
            row_color = "#ffffbf"
        elif 90.0 <= row["Completeness"] < 95.0:
            row_color = "#fdae61"
        else:
            row_color = "#d73027"

        field_str = [f"[{row_color}]{s}" for s in row]
        table.add_row(*field_str)

    if table_output is not None:
        df.to_csv(table_output, sep="\t", header=True, index=False)

    return table


def quality_control_consensus(consensus_results: Path, consensus_subdir: str = "high_freq", table_output: Path = None):

    """ Create a quality control table from the coverage data and consensus sequences """

    coverage_data = {
        sample.name.replace(".txt", ""): sample
        for sample in (consensus_results / "coverage").glob("*.txt")
    }
    fastp_data = {
        sample.name.replace(".json", ""): sample
        for sample in (consensus_results / "quality_control").glob("*.json")
    }

    combined_files = {}
    for assembly in (consensus_results / "consensus").glob("*.consensus.fasta"):
        name = assembly.name.replace(".consensus.fasta", "")
        
        combined_files[name] = SampleFiles(
            assembly=assembly,
            fastp=fastp_data.get(name),
            samtools=coverage_data.get(name)
        )

    samples = []
    for sample, sample_files in combined_files.items():
        all_reads, qc_reads = get_fastp_data(sample_files.fastp)
        aligned_reads, coverage, mean_depth = get_samtools_data(sample_files.samtools)
        completeness, missing = get_consensus_assembly_data(sample_files.assembly)

        qc = SampleQC(
            name=sample,
            reads=all_reads,
            qc_reads=qc_reads,
            aligned_reads=aligned_reads,
            coverage=coverage,
            mean_depth=mean_depth,
            missing_sites=missing,
            completeness=completeness
        )
        samples.append(qc)

    table_freq_title = "".join([s.capitalize() for s in consensus_subdir.name.split("_")])
    table = create_rich_table(samples, title=f"Monkeypox QC ({table_freq_title})", table_output=table_output)

    rprint(table)



@dataclass
class SampleDistance:
    patient: str
    samples: int
    within_median: float


def snp_distance(dist: Path):
    """
    Compute median SNP distance within and between patients
    Sample identifiers conform to Mona's scheme: MPX_A_1 etc.
    """

    dist_mat = pandas.read_csv(dist, index_col=0)

    # Replace column and index names with extracted patient identifier

    patients = [c.split(".")[0].replace("Consensus_", "").split("_")[1] for c in dist_mat.columns]

    dist_mat.index = patients
    dist_mat.columns = patients

    dist_upper = dist_mat.mask(np.triu(np.ones(dist_mat.shape, dtype=np.bool_)))

    melted = pandas.DataFrame(dist_upper).reset_index().melt('index')

    within_patients = melted[melted['index'] == melted['variable']].dropna()
    between_patients = melted[melted['index'] != melted['variable']].dropna()

    # Within patient with only single isolate comparison to itself is already NaN and excluded

    df = pandas.concat([within_patients, between_patients])
    df['comparison'] = ['within' for _ in within_patients.iterrows()] + ['between' for _ in between_patients.iterrows()]\
    
    df.columns = ['patient1', 'patient2', 'distance', 'comparison']

    fig, ax = plt.subplots(
        nrows=1, ncols=1, figsize=(14, 10)
    )
    
    sns.set_style('white')

    p = sns.boxplot(x="distance", y="comparison", data=df, palette="colorblind", linewidth=2.5, ax=ax)
    sns.stripplot(x="distance", y="comparison", data=df, color="darkgray", alpha=0.8, jitter=0.3, size=8, ax=ax)
 

    p.set_xticks(range(int(df['distance'].max())+1))
    p.set_xticklabels(range(int(df['distance'].max())+1))
    plt.xlabel(f"\nSNP distance", fontsize=12, fontweight="bold")
    plt.ylabel(f"Patient isolates\n", fontsize=12, fontweight="bold")
    plt.tight_layout()
    
    fig.savefig("test.png")
