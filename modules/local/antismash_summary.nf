process ANTISMASH_SUMMARY {
    tag "${meta.id}"
    label 'process_low'

    container "microbiome-informatics/mgnify-pipelines-toolkit:1.2.0--htslib"

    input:
    tuple val(meta), path(antismash_gff)

    output:
    tuple val(meta), path("*_summary.tsv.gz")    , emit: antismash_summary
    tuple val(meta), path("*_summary.tsv.gz.gzi"), emit: antismash_summary_index
    path "versions.yml"                          , emit: versions

    script:
    """
    summarise_antismash_bgcs \\
        --antismash-gff ${antismash_gff} \\
        --output ${antismash_gff.simpleName}_summary.tsv

    bgzip -@${task.cpus} --index ${antismash_gff.simpleName}_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """

    stub:
    """
    touch ${antismash_gff.simpleName}_summary.tsv

    bgzip -@${task.cpus} --index ${antismash_gff.simpleName}_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """
}
