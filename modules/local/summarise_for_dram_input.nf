process SUMMARISE_FOR_DRAM_INPUT {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pandas:2.2.3--e136a7b7218cc69c':
        'community.wave.seqera.io/library/pandas:2.2.3--9b034ee33172d809' }"

    input:
    tuple val(meta), path(fasta), path(ko_per_contigs_tsv), path(interproscan_tsv), path(dbcan_overview)

    output:
    tuple val(meta), path("${prefix}_summary_for_dram.tsv"), emit: dram_summary
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # For every assembly analysis, this script extracts:
    # - Columns X and Y from the interpro summary (namely, Pfam ID and description)
    # - A consensus from dbcan_overview.txt for CAZy families
    # - Kegg Orthologs IDs and description
    # It then produces a tsv table for dram distill to generate tabular and visual annotation summaries

    summarise_for_dram.py \\
        --prefix ${prefix} \\
        --fasta ${fasta} \\
        --ko-per-contig ${ko_per_contigs_tsv} \\
        --interproscan ${interproscan_tsv} \\
        --dbcan-overview ${dbcan_overview}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_dram_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
