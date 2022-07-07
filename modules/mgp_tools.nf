process ExtractAligned {

    label "mgp_tools"
    tag { "$id : $idx_name" }

    publishDir "$params.outdir/$params.stage/$params.subdir", mode: "copy", pattern: "${id}_${idx_name}.json"
    publishDir "$params.outdir/$params.stage/$params.subdir/fastq", mode: "symlink", pattern: "${id}_${idx_name}_aligned*.fastq"  // always to fastq sub-subdir

    input:
    tuple val(id), file(forward), file(reverse)
    tuple val(id_idx), val(idx_name), file(alignment)  // can be .txt extension for read identifiers (recognized automatically)

    output:
    tuple val(id), file("${id}_${idx_name}_aligned_1.fastq"), file("${id}_${idx_name}_aligned_2.fastq")
    file("${id}_${idx_name}.json") optional true 

    script:

    if (forward.size() > 0 && reverse.size() > 0) // guard against Needletail error on empty file
        """
        mgp-tools deplete \
            --alignment $alignment \
            --input $forward \
            --input $reverse \
            --output ${id}_${idx_name}_aligned_1.fastq \
            --output ${id}_${idx_name}_aligned_2.fastq \
            --min-len $params.extract_min_qaln_len \
            --min-mapq $params.extract_min_mapq \
            --report ${id}_${idx_name}.json \
            --output-removed
        """
    else
        """
        touch ${id}_${idx_name}_aligned_1.fastq 
        touch ${id}_${idx_name}_aligned_2.fastq
        """

}