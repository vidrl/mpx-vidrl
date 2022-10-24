import groovy.json.JsonOutput


process ArticNanoq {

    // Reimplementation of `artic gather` using `nanoq` for easier
    // parsing and filtering of sequence data without need for 
    // specific sub-directory structure and with support for 
    // compressed files

    tag { "$id" }
    label "artic_gather"

    publishDir "${params.outdir}/quality_control", mode: "copy", pattern: "${id}.nanoq.json"  // publish all outputs as symlinks
    publishDir "${params.outdir}/quality_control", mode: "symlink", pattern: "${id}.read_lengths.txt"  // publish all outputs as symlinks
    publishDir "${params.outdir}/quality_control", mode: "symlink", pattern: "${id}.read_qualities.txt"  // publish all outputs as symlinks

    input:
    tuple val(id), file(fastq_files)

    output:
    tuple val(id), file("${id}.fastq")
    tuple file("${id}.nanoq.json"), file("${id}.read_lengths.txt"), file("${id}.read_qualities.txt")

    script:
    // First concat then filter, so that report has filtered stats properly recorded
    """
    for f in $fastq_files; do
        nanoq -i \$f >> ${id}.fastq
    done 
    nanoq -i ${id}.fastq -s --json --report ${id}.nanoq.json \
        --read-lengths ${id}.read_lengths.txt --read-qualities ${id}.read_qualities.txt \
        --min-len $params.min_length --max-len $params.max_length --min-qual $params.min_quality
    """

}

process ArticMinion {

    tag { "$id" }
    label "artic_minion"

    publishDir "${params.outdir}/consensus/${id}_outputs", mode: "symlink", pattern: "${id}.*"  // publish all outputs as symlinks
    publishDir "${params.outdir}/consensus/", mode: "copy", pattern: "${id}.consensus.fasta"
    publishDir "${params.outdir}/consensus/", mode: "copy", pattern: "${id}.coverage_mask.txt"
    publishDir "${params.outdir}/consensus/", mode: "copy", pattern: "${id}.pass.vcf.gz"

    input:
    tuple val(id), file(fastq)
    tuple file(scheme_base_dir), val(scheme_dir)  // bed file not used just for easier primer scheme passing in main.nf

    output:
    tuple val(id), file("${id}.trimmed.rg.sorted.bam")
    file("${id}.coverage_mask.txt")
    file("${id}.*")

    script:

    """
    artic minion --normalise $params.normalise --threads $task.cpus --read-file $fastq --medaka --medaka-model ${params.medaka_model} --skip-muscle --min-var-depth ${params.medaka_min_depth} --scheme-directory ${scheme_base_dir} ${scheme_dir} ${id}
    """

}

process ArticCoverage {

    tag { "$id" }
    label "artic_covtobed"

    publishDir "${params.outdir}/coverage", mode: "symlink", pattern: "${id}.coverage.*"


    input:
    tuple val(id), file(reg_trimmed_bam)

    output:
    file("${id}.coverage.bed")
    file("${id}.coverage.txt")

    script:

    """
    samtools coverage $reg_trimmed_bam > ${id}.coverage.txt
    covtobed $reg_trimmed_bam > ${id}.coverage.bed
    """

}


process ArticReport {

    publishDir "${params.outdir}/", mode: "copy", pattern: "report.html"

    tag { "Report" }
    label "artic_report"

    input:
    file(coverage_files)
    file(mask_files)
    file(nanoq_files)
    file(scheme_bed)
    file(params_json)

    output:
    file("report.html")

    script:

    if (params.barcodes){
        """
        python $baseDir/modules/artic/report/report.py --scheme $scheme_bed --params $params_json \
            --coverage-glob "*.coverage.bed" \
            --mask-glob "*.coverage_mask.txt" \
            --read-lengths-glob "*.read_lengths.txt" \
            --read-qualities-glob "*.read_qualities.txt" \
            --nanoq-report-glob "*.nanoq.json" \
            --barcodes
        """
    } else {

        """
        python $baseDir/modules/artic/report/report.py --scheme $scheme_bed --params $params_json  \
            --coverage-glob "*.coverage.bed" \
            --mask-glob "*.coverage_mask.txt" \
            --read-lengths-glob "*.read_lengths.txt" \
            --read-qualities-glob "*.read_qualities.txt" \
            --nanoq-report-glob "*.nanoq.json"
        """
    }

}

process ArticParams {

    publishDir "${params.outdir}/", mode: "copy", pattern: "params.json"

    tag { "Parameters" }
    label "artic_report"

    input:
    val(workflow_started)

    output:
    file("params.json")

    script:

    
    def param_data = """
{
    "outdir": "$params.outdir",
    "version": "$params.version",
    "fastq_gather": "$params.fastq_gather",
    "fastq_id": "$params.fastq_id",
    "fastq_dir": "$params.fastq_dir",
    "fastq_ext": "$params.fastq_ext",
    "barcodes": "$params.barcodes",
    "sample_sheet": "$params.sample_sheet",
    "scheme_dir": "$params.scheme_dir",
    "medaka_model": "$params.medaka_model",
    "min_length": "$params.min_length",
    "max_length": "$params.max_length",
    "min_quality": "$params.min_quality",
    "normalise": "$params.normalise",
    "report_title": "$params.report_title",
    "started": "$workflow_started"
}
    """

    """
    cat <<< '
    $param_data
    ' > params.json 
    """
}