process SUMMARISEFORDRAMINPUT {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python:3.13.1--9856f872fdeac74e':
        'community.wave.seqera.io/library/python:3.13.1--d00663700fcc8bcf' }"

    input:
    tuple val(meta), path(root_folder)

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
    # - The best consensus from dbcan_overview.txt for CAZy families
    # - Kegg Orthologs IDs and description
    # It then produces a tsv table for dram distill to generate tabular and visual annotation summaries

    python summarise_for_dram.py -i ${root_folder} -p ${prefix}

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
