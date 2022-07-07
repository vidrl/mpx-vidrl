process Fastp {

    label "fastp"
    tag { id }

    publishDir "$params.outdir/$params.stage/$params.subdir", mode: "copy", pattern: "${id}.json"
    publishDir "$params.outdir/$params.stage/$params.subdir/fastq", mode: "symlink", pattern: "${id}_qc_*.fastq" // always to fastq sub-subdir

    input:
    tuple val(id), file(forward), file(reverse)

    output:
    tuple val(id), file("${id}_qc_1.fastq"), file("${id}_qc_2.fastq")
    file("${id}.json")

    """
    fastp -i $forward -I $reverse -o ${id}_qc_1.fastq -O ${id}_qc_2.fastq --thread $task.cpus --length_required 50 --cut_tail --cut_tail_mean_quality 20 --json ${id}.json
    """

}