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
include { Ivar as IvarHighFrequency }  from './modules/ivar' addParams(
    stage: "consensus",
    subdir: "high_freq",
    ivar_min_freq: params.ivar_min_freq_high
)
include { Ivar as IvarLowFrequency } from './modules/ivar' addParams(
    stage: "consensus",
    subdir: "low_freq",
    ivar_min_freq: params.ivar_min_freq_low
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

// Read quality control, coverage data, consensus assembly, variant calls
workflow qc_variants_assembly {
    take:
        reads
        reference
        gff
    main:
        qc_reads = Fastp(reads)
        aligned_reads = MinimapAlignSortedBam(qc_reads[0], reference)
        coverage = Coverage(aligned_reads)
        if (params.ivar_freq_low) {
            ivar__low = IvarLowFrequency(aligned_reads, reference, gff)
        }
        if (params.ivar_freq_high) {
            ivar_high = IvarHighFrequency(aligned_reads, reference, gff)
        }
}

workflow {

    reads = channel.fromFilePairs(params.fastq, flat: true)
    reference = check_file(params.reference)
    gff = check_file(params.gff)
    qc_variants_assembly(reads, reference, gff)

}

