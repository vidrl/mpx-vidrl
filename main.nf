#!/usr/bin/env nextflow


/* 
vim: syntax=groovy
-*- mode: groovy;-*-

Monkeypox assembly workflow

@authors: Eike Steinig, Mona Taouk

*/

nextflow.enable.dsl=2

include { check_file } from './modules/utils'

include { Fastp } from './modules/fastp' addParams(
    stage: "quality_control",
    subdir: ""
)
include { MinimapAlignCigarPAF as MinimapReferenceAlignment } from './modules/minimap2' addParams(
    stage: "reference_assembly",
    subdir: "alignments"
    align_label: "minimap2"
)
include { ExtractAligned } from './modules/mgp_tools' addParams(
    stage: "reference_assembly",
    subdir: "extractions",
    extract_min_len: params.extract_min_qaln_len,
    extract_min_mapq: params.extract_min_mapq
)
include { Spades as ReferenceAssembly } from './modules/spades' addParams(
    stage: "reference_assembly",
    subdir: "assembly"
)

// Aligns reads to outbreak reference and assembles aligned reads
workflow reference_assembly {
    take:
        reads
        reference
    main:
        aligned_reads = MinimapReferenceAlignment(reads, reference)
        extracted_reads = ExtractAligned(aligned_reads[0], aligned_reads[1])
        assembly = ReferenceAssembly(extracted_reads[0])
    emit:
        assembly[0]
        assembly[1]
}

include { MinimapAlignCigarPAF as MinimapHostAlignment } from './modules/minimap2' addParams(
    stage: "denovo_assembly",
    subdir: "host_alignments"
    align_label: "minimap2_host"
)
include { DepleteAligned } from './modules/mgp_tools' addParams(
    stage: "denovo_assembly",
    subdir: "host_depletion",
    extract_min_len: params.host_deplete_min_qaln_len,
    extract_min_mapq: params.host_deplete_min_mapq
)
include { Spades as DenovoAssembly } from './modules/spades' addParams(
    stage: "denovo_assembly",
    subdir: "assembly"
)

// Depletes human reads and assembles remaining reads
workflow denovo_assembly {
    take:
        reads
        host_index
    main:
        assembly = DenovoAssembly(depleted_reads[0]) 
    emit:
        assembly[0]
        assembly[1]
}

include { MinimapAlignSortedBam } from './modules/minimap2' addParams(
    stage: "consensus_assembly",
    subdir: "alignments"
)

include { IvarConsensus } from './modules/ivar' addParams(
    stage: "consensus_assembly",
    subdir: "consensus"
)

workflow consensus_assembly {
    take:
        reads
        reference
    main:
        aligned_reads = MinimapAlignSortedBam(reads, reference)
        // NB for Mona: We are using only a single reference ever, so it's fine to use
        // it in both processes - if using a channel of multiple references, this might
        // generate wrong combinations of Minimap and Ivar reference files, and it would
        // be better to pass the reference in the alignment output channels directly to 
        // ensure the downstream process gets the right reference (or do some channel
        // operator magic, but it can get messy quickly)
        consensus_assembly = IvarConsensus(aligned_reads, reference)
    emit:
        consensus_assembly
}

workflow {

    reads = channel.fromFilePairs(params.fastq, flat: true)
    
    // All subworkflows use quality controlled reads
    qc_reads = Fastp(reads)

    if (params.host_depletion){
        host_aligned_reads = MinimapHostAlignment(reads, host_index)
        depleted_reads = DepleteAligned(host_aligned_reads[0], host_aligned_reads[1])
        reads = depleted_reads[0]
    } else {
        reads = qc_reads[0]
    }

    if (params.reference_assembly){
        // Reference alignment and assembly of aligned reads
        reference = check_file(params.reference)
        reference_assembly(reads, reference)
    } else if (params.consensus_assembly) {
        // Generate a consensus assembly against the current outbreak reference
        reference = check_file(params.reference)
        consensus_assembly(reads, reference)
    } else if (params.denovo_assembly) {
        // Generate a denovo assembly
        host_index = check_file(params.host_index)
        denovo_assembly(reads, host_index)
    } else {
        error "\nRequired argument mising (--reference_assembly | --consensus_assembly | --denovo_assembly)  [set to `true`]"
    }
}

