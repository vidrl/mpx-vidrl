process SpadesIsolate {

    tag { "$id : $idx_name" }
    label "spades"

    publishDir "$params.outdir/$params.stage/$params.subdir", mode: "symlink", pattern: "${id}.*.fasta"

    input:
    tuple val(id), file(forward), file(reverse)

    output:
    tuple val(id), file(forward), file(reverse)
    tuple val(id), val("${id}.contigs.fasta"), val("${id}.scaffolds.fasta")

    script:

    mem = task.memory.replaceAll("[^0-9]", "")

    """
    spades.py --isolate -1 $forward -2 $reverse -t $task.cpus -m $mem -o working
    cp working/contigs.fasta ${id}.contigs.fasta
    cp working/scaffolds.fasta ${id}.scaffolds.fasta
    """

}

