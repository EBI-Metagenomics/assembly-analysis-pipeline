process SEQKIT_SPLIT2 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.8.1--h9ee0642_0' :
        'biocontainers/seqkit:2.8.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(assembly)
    val(length)
    val(size)

    output:
    tuple val(meta), path("**/*.gz"), emit: assembly
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Validate
    if (length != null && size != null) {
        error("Cannot use both --by-length (${length}) and --by-size (${size}) simultaneously!")
    }
    if (length == null && size == null) {
        error("Neither length nor nparts specified for ${meta.id}. Using seqkit default behavior.")
    }

    def chunk_by_length = (length != null) ? "--by-length ${length}" : ""
    def chunk_in_size   = (size != null) ? "--by-size ${size}" : ""
    """
    seqkit \\
        split2 \\
        $args ${chunk_by_length} ${chunk_in_size} \\
        --threads $task.cpus \\
        ${assembly} \\
        --out-dir ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(echo \$(seqkit 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
