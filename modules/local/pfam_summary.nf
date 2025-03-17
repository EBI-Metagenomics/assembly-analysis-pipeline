process PFAM_SUMMARY {

    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/csvtk:0.31.0--h9ee0642_0':
        'biocontainers/csvtk:0.31.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(interproscan_tsv)

    output:
    tuple val(meta), path("${prefix}_pfam_summary.tsv.gz"), emit: pfam_summary
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # This csvtk minipipeline extracts Pfam rows:
    # . Fix the quotation, sometimes the TSV has quotes where it shouldn't
    # . Picks the Pfam rows.
    # . Extracts the 5th (Pfam accession) and 6th (description) columns.
    # . Counts the frequency of the pfams (we use pfam accession and description because we need the description too)
    # . Adds headers ('pfam', 'description' and 'count') the TSV.
    # . Inverts the columns - we need count to the be first

    csvtk fix-quotes --tabs --no-header-row ${interproscan_tsv} | \\
    csvtk filter2 --tabs --no-header-row --filter '\$4 == "Pfam"' | \\
    csvtk cut --tabs --no-header-row --fields 5,6 | \\
    csvtk freq --tabs --no-header-row --fields 1,2 --reverse --sort-by-freq | \\
    csvtk add-header --tabs --no-header-row --names pfam,description,count | \\
    csvtk cut --tabs --fields count,pfam,description --out-file ${prefix}_pfam_summary.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_pfam_summary.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
    END_VERSIONS
    """
}
