process Spades {

    tag { "$id : $idx_name" }
    label "spades"

    publishDir "$params.outdir/$params.stage/$params.subdir", mode: "symlink", pattern: "*.fasta"

    input:
    tuple val(id), file(forward), file(reverse)

    output:
    tuple val(id), file(forward), file(reverse)
    tuple val(id), val("${id}.contigs.fasta"), val("${id}.scaffolds.fasta")

    script:

    mem = task.memory.toString().replaceAll("[^0-9]", "")

    if (forward.size() > 0 && reverse.size() > 0) // guard against empty file
        """
        spades.py $params.spades_opts -1 $forward -2 $reverse -t $task.cpus -m $mem -o working
        cp working/contigs.fasta ${id}.contigs.fasta
        cp working/scaffolds.fasta ${id}.scaffolds.fasta
        """
    else
        """
        touch ${id}.contigs.fasta
        touch ${id}.scaffolds.fasta
        """
}

