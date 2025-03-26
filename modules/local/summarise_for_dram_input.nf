process SUMMARISE_FOR_DRAM_INPUT {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'depot.galaxyproject.org/singularity/pandas-2.2.1':
        'biocontainers/pandas-2.2.1' }"

    input:
    tuple val(meta), path(ko_summaries)
    tuple val(meta), path(ko_per_contigs)
    tuple val(meta), path(interpro_summaries)
    tuple val(meta), path(dbcan_overviews)

    output:
    tuple val(meta), path("${prefix}_dram_summary.tsv.gz"), emit: dram_summary
    path "versions.yml"                                   , emit: versions

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

    summarise_for_DRAM.py \\
    --prefix ${prefix} \\
    --ko_summaries ${ko_summaries} \\
    --ko_per_contigs ${ko_per_contigs} \\
    --interpro_summaries ${interpro_summaries} \\
    --dbcan_overviews ${dbcan_overviews}

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
