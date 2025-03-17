process INTERPRO_SUMMARY {

    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/csvtk:0.31.0--h9ee0642_0':
        'biocontainers/csvtk:0.31.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(interproscan_tsv)

    output:
    tuple val(meta), path("${prefix}_intepro_summary.tsv.gz"), emit: interpro_summary
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # This csvtk minipipeline extracts InterProScan - InterPro accession counts:
    # . Fix the quotation, sometimes the TSV has quotes where it shouldn't
    # . Picks the rows with an InterPro hit which are those where the InterPro column is not empty or with "-" character
    # . Extracts the 12th InterPro annotations - accession (e.g. IPR002093) and 13th InterPro annotations - description (e.g. BRCA2 repeat)
    # . Counts the frequency of the IPS accessions (we use accession and description because we need the description too)
    # . Adds headers ('interpro_accession', 'description' and 'count') the TSV.
    # . Inverts the columns - we need count to the be first

    csvtk fix-quotes --tabs --no-header-row ${interproscan_tsv} | \\
    csvtk filter2 --tabs --no-header-row --filter '\$12 != "" && \$12 != "-"' | \\
    csvtk cut --tabs --no-header-row --fields 12,13 | \\
    csvtk freq --tabs --no-header-row --fields 1,2 --reverse --sort-by-freq | \\
    csvtk add-header --tabs --no-header-row --names interpro_accession,description,count | \\
    csvtk cut --tabs --fields count,interpro_accession,description > ${prefix}_intepro_summary.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_intepro_summary.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
    END_VERSIONS
    """
}
