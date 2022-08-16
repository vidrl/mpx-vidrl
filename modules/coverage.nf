process Coverage {

    tag { "$id : $idx_name" }
    label "coverage"

    publishDir "$params.outdir/$params.stage/$params.subdir", mode: "copy", pattern: "${id}.coverage.*"

    input:
    tuple val(id), val(idx_name), file(bam)

    output:
    tuple val(id), val("${id}.coverage.txt"), val("${id}.coverage.bed")

    script:

    """
    samtools coverage $bam > ${id}.coverage.txt
    covtobed $bam > ${id}.coverage.bed
    """

}