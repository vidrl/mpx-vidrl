#!/usr/bin/env nextflow


/* 
vim: syntax=groovy
-*- mode: groovy;-*-

Monkeypox assembly workflow

@authors: Mona Taouk, Eike Steinig

*/

nextflow.enable.dsl=2


// Basic workflow using alignment against Monkeypox reference
// and subsequent assembly of mapped reads

include { check_file } from './modules/utils'

include { MinimapAlignCigarPAF as MinimapReferenceAlignment } from './modules/minimap2' addParams(
    stage: "reference_alignment",
    subdir: "alignments"
)
include { ExtractAligned } from './modules/mgp_tools' addParams(
    stage: "reference_alignment",
    subdir: "extractions",
    extract_min_len: params.extract_min_qaln_len,
    extract_min_mapq: params.extract_min_mapq
)
include { SpadesIsolate } from './modules/spades' addParams(
    stage: "reference_alignment",
    subdir: "assembly"
)

workflow {

    // Staging files and channels

    mpx_ref = check_file(params.reference)
    reads = channel.fromFilePairs(params.fastq, flat: true)

    // Mapping reads to MPX reference and extracting mapped reads

    aligned_reads = MinimapReferenceAlignment(reads, mpx_ref)  // output: id, fwd, rev (input reads) + id, idx_name, alignment (alignment)
    extracted_reads = ExtractAligned(aligned_reads[0], aligned_reads[1])  // output: id, fwd, rev (aligned reads) + report json

    // Assembly of reads aligned against MPX reference (Spades)

    assemblies = SpadesIsolate(extracted_reads[0]) // output: id, fwd, rev (aligned reads) + id, contigs, scaffolds


}

