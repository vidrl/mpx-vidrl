process Coverage {

    tag { "$id : $idx_name" }
    label "coverage"

    publishDir "$params.outdir/$params.stage/$params.subdir", mode: "copy", pattern: "${id}.coverage.*"

    input:
    tuple val(id), val(idx_name), file(bam)

    output:
    tuple val(id), file("${id}.coverage.txt"), file("${id}.coverage.bed")

    script:

    """
    samtools coverage $bam > ${id}.coverage.txt
    covtobed $bam > ${id}.coverage.bed
    """

}