process FILTER_ASSEMBLY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.8.1--h9ee0642_0':
        'biocontainers/seqkit:2.8.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${prefix}_filtered.fasta.gz") , emit: fastx
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix  = task.ext.prefix ?: "${meta.id}"
    """
    seqkit \\
        seq \\
        --min-len ${params.min_contig_length} \\
        --threads $task.cpus \\
        $assembly | \\
    seqkit fx2tab \\
        --threads $task.cpus \\
        --base-content N | \\
    awk '\$3 < 10' | \\
    seqkit tab2fx \\
        --threads $task.cpus \\
        --out-file ${prefix}_filtered.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_filtered.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """
}
