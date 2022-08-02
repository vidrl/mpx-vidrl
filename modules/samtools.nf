process SamtoolsCoverage {

    tag { "$id : $idx_name" }
    label "samtools"

    publishDir "$params.outdir/$params.stage/$params.subdir", mode: "symlink", pattern: "${id}_${idx_name}.txt"

    input:
    tuple val(id), val(idx_name), file(bam)

    script:

    idx_name = index.baseName

    """
    minimap2 -t $task.cpus -ax sr ${index} $forward $reverse > ${id}_${idx_name}.sam
    samtools view -S ${id}_${idx_name}.sam -b | samtools sort - -o ${id}_${idx_name}.sorted.bam
    rm ${id}_${idx_name}.sam
    """

}