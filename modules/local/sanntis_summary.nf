process SANNTIS_SUMMARY {
    tag "${meta.id}"
    label 'process_low'

    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
    //    "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"
    container "community.wave.seqera.io/library/pip_mgnify-pipelines-toolkit:2764846cec7da9cf"

    input:
    tuple val(meta), path(sanntis_gff)

    output:
    tuple val(meta), path("*_summary.tsv"), emit: sanntis_summary
    path "versions.yml"                   , emit: versions

    script:
    """
    summarise_sanntis_bgcs \\
        --sanntis-gff ${sanntis_gff} \\
        --output ${sanntis_gff.simpleName}_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """

    stub:
    """
    touch ${sanntis_gff.simpleName}_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """
}
