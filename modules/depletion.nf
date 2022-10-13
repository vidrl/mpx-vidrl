

process Minimap2HostPaired {

    tag { "$id" }
    label "minimap2_host"

    publishDir "$params.outdir/$params.stage/$params.subdir", mode: "symlink", pattern: "${id}.paf"

    input:
    tuple val(id), file(forward), file(reverse)
    file(index)

    output:
    tuple val(id), file(forward), file(reverse)
    tuple val(id), file("${id}.paf")

    script:

    """
    minimap2 -t $task.cpus -c -x sr ${index} $forward $reverse > ${id}.paf
    """

}

process DepleteHostPaired {

    label "mgp_tools"
    tag { "$id" }

    publishDir "$params.outdir/$params.stage/$params.subdir", mode: "copy", pattern: "${id}.json"
    publishDir "$params.outdir/$params.stage/$params.subdir/fastq", mode: "symlink", pattern: "${id}_depleted*.fastq"  // always to fastq sub-subdir

    input:
    tuple val(id), file(forward), file(reverse)
    tuple val(id_idx), file(alignment)  // can be .txt extension for read identifiers (recognized automatically)

    output:
    tuple val(id), file("${id}_depleted_1.fastq"), file("${id}_depleted_2.fastq")
    file("${id}.json") optional true 

    script:

    if (forward.size() > 0 && reverse.size() > 0) // guard against Needletail error on empty file
        """
        mgp-tools deplete --alignment $alignment --input $forward --input $reverse --output ${id}_depleted_1.fastq --output ${id}_depleted_2.fastq --min-cov $params.deplete_min_cov --min-len $params.deplete_min_len --min-mapq $params.deplete_min_mapq --report ${id}.json
        """
    else
        """
        touch ${id}_depleted_1.fastq 
        touch ${id}_depleted_2.fastq
        """

}


process Minimap2HostSingle {

    tag { "$id : $idx_name" }
    label "minimap2_host"

    publishDir "$params.outdir/$params.stage/$params.subdir", mode: "symlink", pattern: "${id}.paf"

    input:
    tuple val(id), file(reads)
    file(index)

    output:
    tuple val(id), file(reads)
    tuple val(id), file("${id}.paf")

    script:

    """
    minimap2 -t $task.cpus -c -x map-ont ${index} $reads > ${id}.paf
    """

}

process DepleteHostSingle {

    label "mgp_tools"
    tag { "$id : $id_idx" }

    publishDir "$params.outdir/$params.stage/$params.subdir", mode: "copy", pattern: "${id}.json"
    publishDir "$params.outdir/$params.stage/$params.subdir/fastq", mode: "symlink", pattern: "${id}_depleted.fastq"  // always to fastq sub-subdir

    input:
    tuple val(id), file(reads)
    tuple val(id_idx), file(alignment)  // can be .txt extension for read identifiers (recognized automatically)

    output:
    tuple val(id), file("${id}_depleted.fastq")
    file("${id}.json") optional true 

    script:

    if (reads.size() > 0) // guard against Needletail error on empty file
        """
        mgp-tools deplete --alignment $alignment --input $reads --output ${id}_depleted.fastq --min-cov $params.deplete_min_cov --min-len $params.deplete_min_len --min-mapq $params.deplete_min_mapq --report ${id}.json
        """
    else
        """
        touch ${id}_depleted.fastq 
        """

}
