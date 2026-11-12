process FILTER_GFF_BY_CONTIGS {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python:3.13.1--9856f872fdeac74e':
        'community.wave.seqera.io/library/python:3.13.1--d00663700fcc8bcf' }"

    input:
    tuple val(meta), path(contigs_chunk), path(gff)

    output:
    // The GFF is not emitted if no CDS belongs to the contigs of the chunk
    tuple val(meta), path(contigs_chunk), path("${prefix}.gff.gz"), emit: filtered, optional: true
    path "versions.yml",                                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    filter_gff_by_contigs.py \\
        --contigs ${contigs_chunk} \\
        --gff ${gff} \\
        --output ${prefix}.gff.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "##gff-version 3" | gzip > ${prefix}.gff.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
