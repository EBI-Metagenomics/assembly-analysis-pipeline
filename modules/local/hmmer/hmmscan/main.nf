process HMMER_HMMSCAN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.4--hdbdd923_1' :
        'biocontainers/hmmer:3.4--hdbdd923_1' }"

    input:
    tuple val(meta), path(hmmfile), path(seqdb), val(write_target), val(write_domain)

    output:
    tuple val(meta), path('*.tbl.gz')   , emit: target_summary, optional: true
    tuple val(meta), path('*.domtbl.gz'), emit: domain_summary, optional: true
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args   ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    target_summary = write_target    ? "--tblout ${prefix}.tbl" : ''
    domain_summary = write_domain    ? "--domtblout ${prefix}.domtbl" : ''
    """
    HMMDB=`find -L ./ -name "*.hmm"`

    hmmscan \\
        --noali \\
        --cut_ga \\
        $args \\
        --cpu $task.cpus \\
        $target_summary \\
        $domain_summary \\
        \$HMMDB \\
        $seqdb > /dev/null

    gzip --no-name ${write_target ? '*.tbl' : ''} \\
        ${write_domain ? '*.domtbl' : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmscan -h | grep -o '^# HMMER [0-9.]*' | sed 's/^# HMMER *//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ${write_align ? "touch ${prefix}.sto" : ''} \\
    ${write_target ? "touch ${prefix}.tbl" : ''} \\
    ${write_domain ? "touch ${prefix}.domtbl" : ''}

    gzip --no-name  ${write_target ? '*.tbl' : ''} \\
        ${write_domain ? '*.domtbl' : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmscan -h | grep -o '^# HMMER [0-9.]*' | sed 's/^# HMMER *//')
    END_VERSIONS
    """
}
