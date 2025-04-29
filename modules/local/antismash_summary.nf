process ANTISMASH_SUMMARY {
    tag "${meta.id}"
    label 'process_low'

    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
    //    "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    container "community.wave.seqera.io/library/pip_mgnify-pipelines-toolkit:2764846cec7da9cf"

    input:
    tuple val(meta), path(antismash_gff)

    output:
    tuple val(meta), path("*_summary.tsv.gz"), emit: antismash_summary
    path "versions.yml"                      , emit: versions

    script:
    """
    summarise_antismash_bgcs \\
        --antismash-gff ${antismash_gff} \\
        --output ${antismash_gff.simpleName}_summary.tsv

    gzip ${antismash_gff.simpleName}_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """

    stub:
    """
    touch ${antismash_gff.simpleName}_summary.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """
}
