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
=======================
A R T I C - P A R A M S
=======================
*/


include { validate_primer_scheme } from './modules/artic/utils'
include { get_fastq_files as get_fastq_files_artic } from './modules/artic/utils'

include { ArticCovtobed } from './modules/artic/artic'
include { ArticReport } from './modules/artic/artic'
include { ArticNanoq } from './modules/artic/artic' addParams(
    min_length: params.artic.min_length,
    max_length: params.artic.max_length,
    min_quality: params.artic.min_quality
)
include { ArticMinion } from './modules/artic/artic' addParams(
    normalise: params.artic.normalise,
    medaka_model: params.artic.medaka_model
)
include { ArticParams } from './modules/artic/artic' addParams(
    outdir: params.artic.outdir,
    version: params.artic.version,
    fastq_dir: params.artic.fastq_dir,
    fastq_ext: params.artic.fastq_ext,
    sample_sheet: params.artic.sample_sheet,
    scheme_dir: params.artic.scheme_dir,
    medaka_model: params.artic.medaka_model,
    min_length: params.artic.min_length,
    max_length: params.artic.max_length,
    min_quality: params.artic.min_quality,
    normalise: params.artic.normalise,
    report_title: params.artic.report_title,
    started: params.artic.workflow_started
)

workflow mpxv_artic {

    started = String.format('%tF %<tH:%<tM', java.time.LocalDateTime.now())

    println("""
    ====================
    Pipeline parameters
    ====================

    outdir:           $params.outdir
    version           $params.version
    
    ====================
    ARTIC amplicon [ONT]
    ====================

    sample_sheet:     $params.artic.sample_sheet
    scheme_dir:       $params.artic.scheme_dir

    fastq_dir:        $params.artic.fastq_dir
    fastq_ext:        $params.artic.fastq_ext

    medaka_model:     $params.artic.medaka_model
    min_length:       $params.artic.min_length
    max_length:       $params.artic.max_length
    min_quality:      $params.artic.min_quality
    normalise:        $params.artic.normalise
    report_title:     $params.artic.report_title
    """)

    if (!params.artic.medaka_model){
        println("Please provide a Medaka model (--medaka_model)")
        System.exit(1)
    }

    if (!params.artic.scheme_dir){
        println("Please provide a primer scheme directory (--scheme_dir)")
        System.exit(1)
    }

    (primer_scheme, primer_bed) = validate_primer_scheme(params.artic.scheme_dir)

    println("Primer scheme directory: ${primer_scheme[0]} (scheme: ${primer_scheme[1]})")
    
    fastq_files = get_fastq_files_artic(
        null, 
        null, 
        params.artic.fastq_dir, 
        params.artic.fastq_ext, 
        null, 
        params.artic.sample_sheet
    )

    artic_params = ArticParams(
        started
    )
    artic_nanoq = ArticNanoq(
        fastq_files
    )
    artic_medaka = ArticMinion(
        artic_nanoq[0], 
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


include { check_file } from './modules/mpxv/utils'
// include { get_fastq_files as get_fastq_files_twist } from './modules/mpxv/utils'

include { Fastp } from './modules/mpxv/mpxv' addParams(
    stage: "quality_control",
    subdir: ""
)
include { MinimapAlignSortedBam } from './modules/mpxv/mpxv' addParams(
    stage: "alignments",
    subdir: ""
)
include { Ivar as IvarHighFrequency }  from './modules/mpxv/mpxv' addParams(
    stage: "consensus",
    subdir: "high_freq",
    ivar_min_freq: params.ivar_min_freq_high
)
include { Ivar as IvarLowFrequency } from './modules/mpxv/mpxv' addParams(
    stage: "consensus",
    subdir: "low_freq",
    ivar_min_freq: params.ivar_min_freq_low
)
include { Coverage } from './modules/mpxv/mpxv' addParams(
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

workflow mpxv_twist {


    started = String.format('%tF %<tH:%<tM', java.time.LocalDateTime.now())

    println("""
    ====================
    Pipeline parameters
    ====================

    outdir:           $params.outdir
    version           $params.version
    

    ===========================
    TWIST enrichment [Illumina]
    ===========================

    sample_sheet:             $params.twist.sample_sheet
    reference:                $params.twist.reference

    fastq_dir:                $params.twist.fastq_dir
    fastq_ext:                $params.twist.fastq_ext

    ivar_ref_gff:             $params.twist.ivar_ref_gff
    ivar_min_qual:            $params.twist.ivar_min_qual
    ivar_min_depth:           $params.twist.ivar_min_depth
    ivar_fill_char:           $params.twist.ivar_fill_char
    ivar_mpileup_args:        $params.twist.ivar_mpileup_args
    ivar_mpileup_max_depth:   $params.twist.ivar_mpileup_max_depth

    """)

    gff = check_file(params.gff)
    reference = check_file(params.reference)

    reads = channel.fromFilePairs(params.fastq, flat: true)

    qc_variants_assembly(reads, reference, gff)

}

