process CONCATENATE_DBCAN_HMMOUT {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/csvtk:0.31.0--h9ee0642_0' :
        'biocontainers/csvtk:0.31.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(tsv, name: 'inputs/tsv*/*')

    output:
    tuple val(meta), path("${prefix}.tsv.gz"), emit: concatenated_tsv
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # The TSV file contains an extra column that appears to be a repeated header.
    # To address this, we first use the `csvtk FIX` command, which adds an empty column
    # to the existing TSV structure. After that, we utilize the `csvtk cat` command
    # to concatenate the modified TSV files, ensuring the output is correctly formatted.

    csvtk fix --tabs inputs/tsv*/** | \\
    csvtk \\
        concat \\
        --num-cpus $task.cpus \\
        --tabs \\
        --out-file ${prefix}.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
    END_VERSIONS
    """

    stub:
    """
    touch ${prefix}.tsv
    gunzip ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
    END_VERSIONS
    """
}
