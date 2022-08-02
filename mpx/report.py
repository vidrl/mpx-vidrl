"""
Monkeypox assembly report
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional
import json

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
    completeness: Optional[float]

def get_fastp_data(file: Path):

    with file.open() as infile:
        fastp_data = json.loads(infile)

    all_reads = fastp_data["summary"]["before_filtering"]["total_reads"]
    qc_reads = fastp_data["summary"]["after_filtering"]["total_reads"]

    return all_reads // 2, qc_reads // 2  # Illumina PE


def get_samtools_data(file: Path):

    content = file.open().readlines()[1].strip().split("\t")
    return content[3], content[5], content[6]  # numreads, coverage, meandepth

def get_assembly_data(file: Path):
    
    pass

def quality_control_consensus(consensus_results: Path):

    """ Create a quality control table from the coverage data and consensus sequences """

    coverage_data = {sample.name.replace(".json", ""): sample for sample in (consensus_results / "coverage").glob("*.txt")}
    fastp_data = {sample.name.replace(".json", ""): sample for sample in (consensus_results / "quality_control").glob("*.json")}

    combined_files = {}
    for assembly in (consensus_results / "consensus_assemblies" / "consensus").glob("*.consensus.fasta"):
        name = assembly.name.replace(".consensus.fasta", "")
        
        combined_files[name] = SampleFiles(
            assembly=assembly,
            fastp=fastp_data.get(name),
            samtools=coverage_data.get(name)
        )
            

    for sample, sample_files in combined_files.items():
        print(f"Processing quality control data for sample: {sample}")

        all_reads, qc_reads = get_fastp_data(sample_files.fastp)
        aligned_reads, coverage, mean_depth = get_samtools_data(sample_files.samtools)
        completeness = get_assembly_data(sample_files.assembly)

