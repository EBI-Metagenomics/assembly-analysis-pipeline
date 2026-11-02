process CONCATENATE_HMMER_TABLES {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4' :
        'biocontainers/pigz:2.3.4' }"

    input:
    tuple val(meta), path(tables)

    output:
    tuple val(meta), path("${prefix}"), emit: file_out
    path 'versions.yml', emit: versions

    script:
    def files = tables.collect { it.toString() }.sort()
    prefix = task.ext.prefix ?: "${meta.id}.tbl.gz"
    """
    zcat ${files.join(' ')} | \\
        awk '!/^#/ && NF' | \\
        pigz --no-name --stdout --processes ${task.cpus} > ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$(pigz --version 2>&1 | sed 's/pigz //g')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.tbl.gz"
    """
    printf '' | pigz --no-name --stdout > ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$(pigz --version 2>&1 | sed 's/pigz //g')
    END_VERSIONS
    """
}
