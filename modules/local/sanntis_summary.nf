process SANNTIS_SUMMARY {
    tag "${meta.id}"
    label 'process_low'

    container "microbiome-informatics/mgnify-pipelines-toolkit:1.2.0--htslib"

    input:
    tuple val(meta), path(sanntis_gff)

    output:
    tuple val(meta), path("*_summary.tsv.gz"),     emit: sanntis_summary,       optional: true
    tuple val(meta), path("*_summary.tsv.gz.gzi"), emit: sanntis_summary_index, optional: true
    path "versions.yml",                           emit: versions

    script:
    """
    # SanntiS produces GFF with just the #gff header if no BGC were found
    if [ "\$(zcat "${sanntis_gff}" | wc -l)" -gt 1 ]; then
        summarise_sanntis_bgcs \\
            --sanntis-gff ${sanntis_gff} \\
            --output ${sanntis_gff.simpleName}_summary.tsv

        bgzip -@${task.cpus} --index ${sanntis_gff.simpleName}_summary.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """

    stub:
    """
    touch ${sanntis_gff.simpleName}_summary.tsv
    bgzip -@${task.cpus} --index ${sanntis_gff.simpleName}_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """
}
