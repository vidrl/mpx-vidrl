process IvarConsensus {

    tag { "$id : $idx_name" }
    label "ivar"

    publishDir "$params.outdir/$params.stage/$params.subdir", mode: "copy", pattern: "*consensus*"

    input:
    tuple val(id), val(idx_name), file(bam)
    file(reference)

    output:
    tuple val(id), file("${id}.consensus.fasta"), file("${id}.consensus.qual.txt")

    script:

    // No `-k` - always add fill character (N) for bases with less than `params.ivar_consensus_min_depth`

    """
    samtools mpileup -d $params.mpileup_max_depth -A -Q 0 $bam | ivar consensus -p ${id}.consensus \
        -q $params.ivar_consensus_min_base_quality \
        -t $params.ivar_consensus_min_frequency \
        -m $params.ivar_consensus_min_depth \
        -n $params.ivar_consensus_fill_char
    mv ${id}.consensus.fa ${id}.consensus.fasta
    """

}