process ANTISMASH_SUMMARY {
    tag "${meta.id}"
    label 'process_low'

    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
    //    "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:1.0.5--pyhdfd78af_0":
        "biocontainers/mgnify-pipelines-toolkit:1.0.5--pyhdfd78af_0" }"

    input:
    tuple val(meta), path(antismash_gff)

    output:
    tuple val(meta), path("*.summary.tsv"), emit: antismash_summary
    path "versions.yml"                   , emit: versions

    script:
    """
    summarise_antismash_bgcs.py
        --antismash-gff ${antismash_gff} \\
        --output ${antismash_gff.simpleName}_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """

    stub:
    """
    touch ${antismash_gff.simpleName}_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """
}
