
process SUMMARISEGOSLIMS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:0.2.0--pyhdfd78af_0':
        'biocontainers/mgnify-pipelines-toolkit:0.2.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(ips)
    tuple val(meta2), path(gaf)
    path go_obo
    path go_banding

    output:
    tuple val(meta), path("*_summary.tsv"), emit: go_summary
    tuple val(meta), path("*_slim.tsv")   , emit: goslim_summary
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ips_uncompressed = ips.name.replace('.gz', '')
    """
    gunzip -c -d ${ips} > ${ips_uncompressed}
    summarise_goslims.py \\
        -go ${go_obo} \\
        -gb ${go_banding} \\
        -i ${ips_uncompressed} \\
        -gaf ${gaf} \\
        -o ${prefix}_summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_summary.tsv
    touch ${prefix}_summary_slim.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """
}
