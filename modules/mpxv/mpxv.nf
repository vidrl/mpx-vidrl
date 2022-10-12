process MpxvMinimapAlignSortedBam {

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
    minimap2 -t $task.cpus -ax sr ${index} $forward $reverse --sam-hit-only > ${id}_${idx_name}.sam
    samtools view -S ${id}_${idx_name}.sam -b | samtools sort - -o ${id}_${idx_name}.sorted.bam
    rm ${id}_${idx_name}.sam
    """

}

process MpxvIvar {

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
    samtools mpileup $params.samtools_mpileup_args -d $params.samtools_mpileup_max_depth -A -B -Q 0 $bam | ivar consensus -p ${id}.consensus \
        -q $params.ivar_min_qual \
        -t $params.ivar_min_freq \
        -m $params.ivar_min_depth \
        -n $params.ivar_fill_char
    
    samtools mpileup $params.samtools_mpileup_args -d $params.samtools_mpileup_max_depth -A -B -Q 0 $bam | ivar variants -p ${id}.variants \
        -q $params.ivar_min_qual \
        -t $params.ivar_min_freq \
        -m $params.ivar_min_depth \
        -r $reference \
        -g $gff

    mv ${id}.consensus.fa ${id}.consensus.fasta
    """

}

process MpxvFastp {

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

process MpxvCoverage {

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