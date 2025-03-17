
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
    tuple val(meta), path("${prefix}/summary.kegg_pathways.tsv.gz"), emit: kegg_pathways
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # This is a custom script used to adjust the structure of the ko per contig tsv:
    # From:
    # ko      contig_id
    # K07082  MGYA1767_1
    # K09458  MGYA1767_1
    # K02990  MGYA1767_1
    # to ->
    # MGYA1767_1    K07082  K09458  K02990

    aggregate_kos_per_contig.py \\
        -i ${ko_contig_tsv} \\
        -o ${prefix}_aggregated_kos_per_contig.tsv

    give_completeness \\
        -i ${prefix}_aggregated_kos_per_contig.tsv \\
        -o ${prefix}

    mv ${prefix}/summary.kegg_pathways.tsv ${prefix}_summary_kegg_pathways.tsv

    gzip ${prefix}_summary_kegg_pathways.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kegg-pathways-completeness: \$(give_completeness --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_summary_kegg_pathways.tsv

    gzip ${prefix}_summary_kegg_pathways.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kegg-pathways-completeness: \$(give_completeness --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
