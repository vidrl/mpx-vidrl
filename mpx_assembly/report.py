"""
Monkeypox assembly report
"""
import pandas
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List
import json
from pyfastx import Fasta
from rich.table import Table
from rich import print as rprint
import pandas

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
    return int(content[3]), float(content[5]), float(content[6])  # numreads, coverage, meandepth


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


def create_rich_table(samples: List[SampleQC], title: str):

    df = pandas.DataFrame(
        [sample.to_list() for sample in samples],
        columns=["Sample", "Reads", "QC Reads", "Alignments", "Coverage", "Mean Depth", "Missing", "Completeness"]
    )

    df = df.sort_values(["Sample", "Completeness", "Coverage", "Mean Depth"])

    print(df)

    table = Table(title=title)
    for cname in df.columns:
        table.add_column(cname, justify="left", no_wrap=False)
    for _, row in df.iterrows():
        field_str = [str(s) for s in row]
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