#!/usr/bin/env nextflow


/* 
vim: syntax=groovy
-*- mode: groovy;-*-

Monkeypox workflow fir Illumina PE from enriched sequencing:
    
    - Quality control
    - Host depletion [optional]
    - Consensus sequences: 3% and 75%
    - Multiple sequence alignment (MSA)
    - Phylogeny
    - Report: QC, SNP distances, Variants, Phylogeny

@authors: Eike Steinig, Mona Taouk
@date: August 2022

*/

nextflow.enable.dsl=2

include { check_file } from './modules/utils'

include { Fastp } from './modules/fastp' addParams(
    stage: "quality_control",
    subdir: ""
)
include { MinimapAlignSortedBam } from './modules/minimap2' addParams(
    stage: "alignments",
    subdir: ""
)
include { IvarConsensus as IvarConsensusHighFrequency }  from './modules/ivar' addParams(
    stage: "consensus",
    subdir: "high_freq",
    ivar_consensus_min_freq: params.ivar_consensus_min_freq_high
)
include { IvarConsensus as IvarConsensusLowFrequency } from './modules/ivar' addParams(
    stage: "consensus",
    subdir: "low_freq",
    ivar_consensus_min_freq: params.ivar_consensus_min_freq_low
)
include { Coverage } from './modules/coverage' addParams(
    stage: "coverage",
    subdir: ""
)

workflow host_depletion {
    take: 
        qc_reads
    main:
        host_index = check_file(params.host_index)
        host_aligned_reads = MinimapHostAlignment(qc_reads[0], host_index)
        // depleted_reads = DepleteAligned(host_aligned_reads[0], host_aligned_reads[1])
    emit:
        depleted_reads[0]
    
}

// Read quality control, consensus assembly, coverage data
workflow qc_consensus_assembly {
    take:
        reads
        reference
    main:
        qc_reads = Fastp(reads)
        aligned_reads = MinimapAlignSortedBam(qc_reads[0], reference)
        coverage = Coverage(aligned_reads)
        consensus_assembly_high = IvarConsensusHighFrequency(aligned_reads, reference)
        consensus_assembly_low = IvarConsensusLowFrequency(aligned_reads, reference)
    emit:
        consensus_assembly_high
        consensus_assembly_low
        coverage
}

workflow {

    reads = channel.fromFilePairs(params.fastq, flat: true)
    reference = check_file(params.reference)
    qc_consensus_assembly(reads, reference)

}

