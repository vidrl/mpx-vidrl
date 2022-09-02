process Ivar {

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