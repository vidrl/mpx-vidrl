#!/usr/bin/env nextflow


/* 
vim: syntax=groovy
-*- mode: groovy;-*-

Monkeypox workflow:

Twist Illumina PE:

    - Quality control
    - Host depletion [optional]
    - Consensus sequences: 60% and 90% frequency thresholds (intra-patient and inter-patient thresholds) [iVar]
    - Variants: > 3% frequency [iVar]

ARTIC ONT amplicons (ARTIC Fieldinformatics):

    - Quality control
    - ARTIC MinION
    - Custom report

@authors: Eike Steinig, Mona Taouk
@date: August 2022

*/

nextflow.enable.dsl=2


/*

*/

/*
=======================
A R T I C - P A R A M S
=======================
*/


include { validate_primer_scheme } from './modules/utils'
include { get_fastq_files as get_fastq_files_artic } from './modules/utils'


include { DepleteHostSingle} from './modules/depletion' addParams(
    stage: "host_depletion",
    subdir: ""
)
include { Minimap2HostSingle } from './modules/depletion' addParams(
    stage: "host_depletion",
    subdir: ""
)

include { ArticCovtobed } from './modules/artic/artic'
include { ArticReport } from './modules/artic/artic'

include { ArticNanoq } from './modules/artic/artic' addParams(
    min_length: params.min_length,
    max_length: params.max_length,
    min_quality: params.min_quality
)
include { ArticMinion } from './modules/artic/artic' addParams(
    normalise: params.normalise,
    medaka_model: params.medaka_model
)
include { ArticParams } from './modules/artic/artic' addParams(
    outdir: params.outdir,
    version: params.version,
    fastq_gather: null,
    fastq_id: null,
    barcodes: null,
    fastq_dir: params.fastq_dir,
    fastq_ext: params.fastq_ext,
    sample_sheet: params.sample_sheet,
    scheme_dir: params.scheme_dir,
    medaka_model: params.medaka_model,
    min_length: params.min_length,
    max_length: params.max_length,
    min_quality: params.min_quality,
    normalise: params.normalise,
    report_title: params.report_title
)

workflow mpxv_artic {

    started = String.format('%tF %<tH:%<tM', java.time.LocalDateTime.now())

    println("""
    ====================
    Pipeline parameters
    ====================

    outdir:           $params.outdir
    version           $params.version
    
    deplete_host:     $params.deplete_host
    host_index:       $params.host_index

    ====================
    ARTIC amplicon [ONT]
    ====================

    sample_sheet:     $params.sample_sheet
    scheme_dir:       $params.scheme_dir

    fastq_dir:        $params.fastq_dir
    fastq_ext:        $params.fastq_ext

    medaka_model:     $params.medaka_model
    min_length:       $params.min_length
    max_length:       $params.max_length
    min_quality:      $params.min_quality
    normalise:        $params.normalise
    report_title:     $params.report_title
    """)

    if (!params.medaka_model){
        println("Please provide a Medaka model (--medaka_model)")
        System.exit(1)
    }

    if (!params.scheme_dir){
        println("Please provide a primer scheme directory (--scheme_dir)")
        System.exit(1)
    }

    (primer_scheme, primer_bed) = validate_primer_scheme(params.scheme_dir)

    println("Primer scheme directory: ${primer_scheme[0]} (scheme: ${primer_scheme[1]})")
    
    fastq_files = get_fastq_files_artic(
        null, 
        null, 
        params.fastq_dir, 
        params.fastq_ext, 
        null, 
        params.sample_sheet
    )

    artic_params = ArticParams(
        started
    )
    artic_nanoq = ArticNanoq(
        fastq_files
    )


    if (params.deplete_host) {
        host_index = check_file(params.host_index)
        host_aligned_reads = Minimap2HostSingle(
            artic_nanoq[0], host_index
        )
        depleted = DepleteHostSingle(
            host_aligned_reads[0], host_aligned_reads[1]
        )
        reads = depleted[0]
    } else {
        reads = artic_nanoq[0]
    }

    artic_medaka = ArticMinion(
        reads, 
        primer_scheme
    )
    artic_coverage = ArticCovtobed(
        artic_medaka[0]
    )

    artic_report = ArticReport(
        artic_coverage | collect,
        artic_medaka[1] | collect,
        artic_nanoq[1] | collect,
        primer_bed,
        artic_params
    )

}

/*
=======================
T W I S T - P A R A M S
=======================
*/


include { check_file } from './modules/utils'
include { get_samples_paired } from './modules/utils'


include { DepleteHostPaired } from './modules/depletion' addParams(
    stage: "host_depletion",
    subdir: ""
)
include { Minimap2HostPaired } from './modules/depletion' addParams(
    stage: "host_depletion",
    subdir: ""
)


include { MpxvFastp } from './modules/mpxv/mpxv'
include { MpxvMinimapAlignSortedBam } from './modules/mpxv/mpxv'

include { MpxvIvar as MpxvIvar90 }  from './modules/mpxv/mpxv' addParams(
    stage: "consensus",
    subdir: "90",
    ivar_min_freq: 0.9
)
include { MpxvIvar as MpxvIvar75 } from './modules/mpxv/mpxv' addParams(
    stage: "consensus",
    subdir: "75",
    ivar_min_freq: 0.75
)
include { MpxvIvar as MpxvIvar60 } from './modules/mpxv/mpxv' addParams(
    stage: "consensus",
    subdir: "60",
    ivar_min_freq: 0.60
)
include { MpxvIvar as MpxvIvar3 } from './modules/mpxv/mpxv' addParams(
    stage: "consensus",
    subdir: "3",
    ivar_min_freq: 0.03
)
include { MpxvIvar as MpxvIvar0 } from './modules/mpxv/mpxv' addParams(
    stage: "consensus",
    subdir: "0",
    ivar_min_freq: 0
)
include { MpxvCoverage } from './modules/mpxv/mpxv' addParams(
    stage: "coverage",
    subdir: ""
)

workflow mpxv_twist {


    started = String.format('%tF %<tH:%<tM', java.time.LocalDateTime.now())

    println("""
    ====================
    Pipeline parameters
    ====================

    outdir:           $params.outdir
    version           $params.version
    
    deplete_host:     $params.deplete_host
    host_index:       $params.host_index

    ===========================
    TWIST enrichment [Illumina]
    ===========================

    sample_sheet:             $params.sample_sheet
    fastq_dir:                $params.fastq_dir

    ivar_ref_gff:             $params.ivar_ref_gff
    ivar_min_qual:            $params.ivar_min_qual
    ivar_min_depth:           $params.ivar_min_depth
    ivar_fill_char:           $params.ivar_fill_char
    ivar_mpileup_args:        $params.ivar_mpileup_args
    ivar_mpileup_max_depth:   $params.ivar_mpileup_max_depth

    """)

    gff = check_file(params.gff)
    reference = check_file(params.reference)

    reads = get_samples(params.fastq_dir, params.sample_sheet)
    
    reference = check_file(params.reference)
    gff = check_file(params.gff)

    reads = MpxvFastp(reads)

    if (params.deplete_host) {
        host_index = check_file(params.host_index)
        host_aligned_reads = Minimap2HostPaired(
            reads[0], host_index
        )
        reads = DepleteHostPaired(
            host_aligned_reads[0], host_aligned_reads[1]
        )
    }
    
    aligned_reads = MpxvMinimapAlignSortedBam(reads[0], reference)
    coverage = MpxvCoverage(aligned_reads)

    MpxvIvar90(aligned_reads, reference, gff)
    MpxvIvar75(aligned_reads, reference, gff)
    MpxvIvar60(aligned_reads, reference, gff)
    MpxvIvar3(aligned_reads, reference, gff)
    MpxvIvar0(aligned_reads, reference, gff)
}

