
process KEGGPATHWAYSCOMPLETENESS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kegg-pathways-completeness:1.3.0--pyhdfd78af_0':
        'biocontainers/kegg-pathways-completeness:1.3.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(ko_contig_tsv)

    output:
    tuple val(meta), path("${prefix}/summary.kegg_pathways.tsv"), emit: kegg_pathways
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // Handle compressde files
    """
    gunzip -c ${ko_contig_tsv} > ${prefix}_ko_contig.tsv

    give_completeness \\
        -i ${prefix}_ko_contig.tsv \\
        -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kegg-pathways-completeness: \$(give_completeness --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.kegg_contigs.tsv
    touch ${prefix}.kegg_pathways.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kegg-pathways-completeness: \$(give_completeness --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
