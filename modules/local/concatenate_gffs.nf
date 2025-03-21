process CONCATENATE_GFFS {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python:3.13.1--9856f872fdeac74e':
        'community.wave.seqera.io/library/python:3.13.1--d00663700fcc8bcf' }"

    input:
    tuple val(meta), path(gffs, name: "gffs/?.gff.gz")

    output:
    tuple val(meta), path('*concatenated.gff') , emit: concatenated_gff
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    concatenate_gffs.py \\
    --gffs gffs/*.gff.gz  \\
    --output ${prefix}_concatenated.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '##gff-version 3' > ${prefix}_concatenated.gff
    gzip ${prefix}_concatenated.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
