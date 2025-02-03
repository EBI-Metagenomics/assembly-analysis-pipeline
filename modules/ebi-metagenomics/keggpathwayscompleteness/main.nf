
process KEGGPATHWAYSCOMPLETENESS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kegg-pathways-completeness:1.3.0--pyhdfd78af_0':
        'biocontainers/kegg-pathways-completeness:1.3.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(proteins_fasta), path(hmmscan_domtbl)

    output:
    // TODO: this contig file is not generated for the example
    // tuple val(meta), path("*_contigs.tsv") , emit: kegg_contigs
    tuple val(meta), path("${prefix}/summary.kegg_pathways.tsv"), emit: kegg_pathways
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // Handle compressde files
    def is_proteins_fasta_compressed = proteins_fasta.name.endsWith(".gz")
    def is_hmmscan_domtbl_compresssed = hmmscan_domtbl.name.endsWith(".gz")
    def proteins_fasta_filename = proteins_fasta.name.replace(".gz", "")
    def hmmscan_domtbl_filename = hmmscan_domtbl.name.replace(".gz", "")
    """
    if [ "${is_proteins_fasta_compressed}" == "true" ]; then
        gzip -c -d ${proteins_fasta} > ${proteins_fasta_filename}
    fi
    if [ "${is_hmmscan_domtbl_compresssed}" == "true" ]; then
        gzip -c -d ${hmmscan_domtbl} > ${hmmscan_domtbl_filename}
    fi

    make_hmmer_table_tab_separated \\
        -i ${hmmscan_domtbl_filename} \\
        -o ${prefix}_kos.tsv

    parse_hmmer_tab_separated_table -t hmmscan \\
        -i ${prefix}_kos.tsv \\
        -f ${proteins_fasta_filename} \\
        -o ${prefix}_per_contig_kos

    give_completeness \\
        -i ${prefix}_per_contig_kos/${prefix}_kos.tsv_parsed \\
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
