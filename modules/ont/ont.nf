process OntMinimapAlignSortedBam {

    tag { "$id : $idx_name" }
    label "minimap2"

    publishDir "$params.outdir/$params.stage/$params.subdir", mode: "symlink", pattern: "${id}_${idx_name}.sorted.bam"

    input:
    tuple val(id), file(reads)
    file(index)

    output:
    tuple val(id), val(idx_name), file("${id}_${idx_name}.sorted.bam")

    script:

    idx_name = index.baseName

    """
    minimap2 -t $task.cpus -ax map-ont ${index} $reads --sam-hit-only > ${id}_${idx_name}.sam
    samtools view -S ${id}_${idx_name}.sam -b | samtools sort - -o ${id}_${idx_name}.sorted.bam
    rm ${id}_${idx_name}.sam
    """

}

process OntPrimerTrim {

    tag { "$id : $idx_name" }
    label "ivar"


    publishDir "$params.outdir/$params.stage/$params.subdir", mode: "symlink", pattern: "${id}_${idx_name}.trimmed.bam"

    input:
    tuple val(id), val(idx_name), file(bam)
    file(primer_bed)

    output:
    tuple val(id), val(idx_name), file("${id}_${idx_name}.trimmed.bam")

    """
    ivar trim -b $primer_bed -p ${id}_${idx_name}.trimmed -i $bam -q 10 -m 100 -s 4
    """

}

process OntIvar {

    tag { "$id : $idx_name" }
    label "ivar"

    publishDir "$params.outdir/$params.stage/$params.subdir", mode: "copy", pattern: "*.consensus.fasta"
    publishDir "$params.outdir/$params.stage/$params.subdir", mode: "copy", pattern: "*.variants.tsv"

    input:
    tuple val(id), val(idx_name), file(bam)
    file(reference)
    file(gff)

    output:
    tuple val(id), file("${id}.consensus.fasta"), file("${id}.variants.tsv")

    script:

    """
    samtools mpileup $params.ivar_mpileup_args -d $params.ivar_mpileup_max_depth -A -B -Q 0 $bam | ivar consensus -p ${id}.consensus \
        -q $params.ivar_min_qual \
        -t $params.ivar_min_freq \
        -m $params.ivar_min_depth \
        -n $params.ivar_fill_char
    
    samtools mpileup $params.ivar_mpileup_args -d $params.ivar_mpileup_max_depth -A -B -Q 0 $bam | ivar variants -p ${id}.variants \
        -q $params.ivar_min_qual \
        -t $params.ivar_min_freq \
        -m $params.ivar_min_depth \
        -r $reference \
        -g $gff

    mv ${id}.consensus.fa ${id}.consensus.fasta
    """

}


process OntCoverage {

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