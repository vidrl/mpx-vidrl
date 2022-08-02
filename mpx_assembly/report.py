"""
Monkeypox assembly report
"""

import json
import numpy

import pandas
from pathlib import Path
from pyfastx import Fasta
from rich.table import Table
from rich import print as rprint
from dataclasses import dataclass
from typing import Optional, List
import numpy as np
from statistics import median

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


def get_fastp_data(file: Path) -> (int, int):
    """
    Get fastp data - divive by two for paired-end reads
    """
    with file.open() as infile:
        fastp_data = json.load(infile)

    all_reads = fastp_data["summary"]["before_filtering"]["total_reads"]
    qc_reads = fastp_data["summary"]["after_filtering"]["total_reads"]

    return all_reads // 2, qc_reads // 2  # Illumina PE


def get_samtools_data(file: Path) -> (int, float, float):
    """
    Get samtools coverage data
    """
    content = file.open().readlines()[1].strip().split("\t")
    return int(content[3]), round(float(content[5]), 4), round(float(content[6]), 6)  # numreads, coverage, meandepth


def get_consensus_assembly_data(file: Path) -> (float or None, int):
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


def create_rich_table(samples: List[SampleQC], title: str, patient_id: bool = True):

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
        rprint(sorted_samples)
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
            row_color = "green3 dim"
        elif 99.0 <= row["Completeness"] < 99.9:
            row_color = "gold3 dim"
        elif 95.0 <= row["Completeness"] < 99.0:
            row_color = "orange3 dim"
        else:
            row_color = "red3 dim"

        field_str = [f"[{row_color}]{s}" for s in row]
        table.add_row(*field_str)
    return table


def quality_control_consensus(consensus_results: Path):

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
    for assembly in (consensus_results / "consensus_assembly" / "consensus").glob("*.consensus.fasta"):
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

    table = create_rich_table(samples, title="Monkeypox QC")
    rprint(table)


@dataclass
class SampleDistance:
    patient: str
    samples: int
    within_median: float


def snp_distance(dist: Path):
    """
    Compute median SNP distance within and between patients
    SAmple identifiers conform to Mona's scheme: MPX_A_1 etc.
    """

    dist_mat = pandas.read_csv(dist, index_col=0)

    # Replace column and index names with extracted patient identifier

    patients = [c.split(".")[0].replace("Consensus_", "").split("_")[1] for c in dist_mat.columns]

    dist_mat.index = patients
    dist_mat.columns = patients

    dist_lower = dist_mat.mask(np.triu(np.ones(dist_mat.shape, dtype=np.bool_)))

    patients_unique = sorted(list(set(patients)))

    between_data = []
    within_data = []
    for patient in patients_unique:

        # Within patient distances
        within_patient = dist_lower.loc[patient, patient]
        if isinstance(within_patient, numpy.float64) and numpy.isnan(within_patient):
            within_median = -1
        else:
            distances = [v for v in within_patient.values.flatten() if not np.isnan(v)]
            within_median = median(distances)

        within_data.append([patient, patient, within_median])

        rprint(
            f"Within patient [red]{patient}[/red] (n = {patients.count(patient)}) "
            f"median SNP distance: [yellow]{within_median}[/yellow]"
        )

        # Between patient distances
        other_patients = [p for p in patients_unique if p != patient]

        for other_patient in other_patients:
            between_patients = dist_lower.loc[patient, other_patient]
            rprint(f"[red]{patient} <--> {other_patient}[/red]")
            # Ignore if all nan, the other combination will have the values:
            if isinstance(between_patients, numpy.float64) and numpy.isnan(between_patients):
                # Single isolate vs. single isolate where value is nan
                median_between = np.nan
            elif isinstance(between_patients, numpy.float64):
                # Single isolate vs. single isolate where value is present
                median_between = median([between_patients])
            else:
                nan_check = np.isnan(between_patients.values).all()
                if nan_check:
                    median_between = np.nan
                else:
                    median_between = median([v for v in between_patients.values.flatten()])

            rprint(f"Between patient median SNP distance: [yellow]{median_between}[/yellow]")
            if not np.isnan(median_between):
                # Sort the patient identifiers (sortable) to fill only single trriangle of matrix
                combo = sorted([patient, other_patient]) + [median_between]
                between_data.append(combo)

    long_form_upper = sorted(between_data, key=lambda x: x[0])
    df = pandas.DataFrame(np.nan, index=patients_unique, columns=patients_unique)

    # Fill in upper triangle with between patient distances
    for data in long_form_upper:
        df.loc[data[0], data[1]] = data[2]

    # Fill in diagonal with within patient distances
    for data in within_data:
        df.loc[data[0], data[1]] = data[2]

    print(df)