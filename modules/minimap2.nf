process MinimapAlignCigarPAF {

    tag { "$id : $idx_name" }
    label "$params.align_label"

    publishDir "$params.outdir/$params.stage/$params.subdir", mode: "symlink", pattern: "${id}_${idx_name}.paf"

    input:
    tuple val(id), file(forward), file(reverse)
    file(index)

    output:
    tuple val(id), file(forward), file(reverse)
    tuple val(id), val(idx_name), file("${id}_${idx_name}.paf")

    script:

    idx_name = index.baseName

    """
    minimap2 -t $task.cpus -c -x sr ${index} $forward $reverse > ${id}_${idx_name}.paf
    """

}

process MinimapAlignSortedBam {

    tag { "$id : $idx_name" }
    label "minimap2"

    publishDir "$params.outdir/$params.stage/$params.subdir", mode: "symlink", pattern: "${id}_${idx_name}.sorted.bam"

    input:
    tuple val(id), file(forward), file(reverse)
    file(index)

    output:
    tuple val(id), val(idx_name), file("${id}_${idx_name}.sorted.bam")

    script:

    idx_name = index.baseName

    """
    minimap2 -t $task.cpus -ax sr ${index} $forward $reverse > ${id}_${idx_name}.sam
    samtools view -S ${id}_${idx_name}.sam -b | samtools sort - -o ${id}_${idx_name}.sorted.bam
    rm ${id}_${idx_name}.sam
    """

}