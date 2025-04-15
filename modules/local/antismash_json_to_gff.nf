process ANTISMASH_JSON_TO_GFF {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}':
        'biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}' }"

    input:
    tuple val(meta), path(antismash_json)

    output:
    tuple val(meta), path("*.gff.gz"), emit: antismash_gff
    path "versions.yml"              , emit: versions

    script:
    """
    # TODO: fix this in the toolkit
    python /usr/local/lib/python3.11/site-packages/mgnify_pipelines_toolkit/analysis/assembly/antismash_gff_builder.py \\
        --input ${antismash_json} \\
        --output ${antismash_json.simpleName}.gff

    gzip ${antismash_json.simpleName}.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """

    stub:
    """
    touch ${antismash_json.simpleName}.gff
    gzip ${antismash_json.simpleName}.gff
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """
}
