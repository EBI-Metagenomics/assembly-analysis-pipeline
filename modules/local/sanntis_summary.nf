process SANNTIS_SUMMARY {
    tag "${meta.id}"
    label 'process_low'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}"
        : "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}"}"

    input:
    tuple val(meta), path(sanntis_gff)

    output:
    tuple val(meta), path("*_summary.tsv.gz"), emit: sanntis_summary, optional: true
    path "versions.yml",                       emit: versions

    script:
    """
    # SanntiS produces GFF with just the #gff header if no BGC were found
    if [ "\$(zcat "${sanntis_gff}" | wc -l)" -gt 1 ]; then
        summarise_sanntis_bgcs \\
            --sanntis-gff ${sanntis_gff} \\
            --output ${sanntis_gff.simpleName}_summary.tsv

        gzip ${sanntis_gff.simpleName}_summary.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """

    stub:
    """
    touch ${sanntis_gff.simpleName}_summary.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """
}
