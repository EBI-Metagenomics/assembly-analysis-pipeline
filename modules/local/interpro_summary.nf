process INTERPRO_SUMMARY {

    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/csvtk:0.31.0--h9ee0642_0':
        'biocontainers/csvtk:0.31.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(interproscan_tsv)

    output:
    tuple val(meta), path("${prefix}_intepro_summary.tsv"), emit: interpro_summary
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # This csvtk minipipeline extracts InterProScan - IPS counts:
    # 1. Picks the rows with an IPS hit
    # 2. Extracts the 12th InterPro annotations - accession (e.g. IPR002093) and 13th InterPro annotations - description (e.g. BRCA2 repeat)xq
    # 3. Counts the frequency of the IPS accessions (we use accession and description because we need the description too)
    # 4. Adds headers ('interpro_accession', 'description' and 'count') the TSV.
    # 5. Inverts the columns - we need count to the be first

    csvtk filter2 --tabs --no-header-row --filter '\$12 != ""' ${interproscan_tsv} | \\
    csvtk cut --tabs --no-header-row --fields 12,13 | \\
    csvtk freq --tabs --no-header-row --fields 1,2 -n | \\
    csvtk add-header --tabs --no-header-row --names interpro_accession,description,count | \\
    csvtk cut --tabs --fields count,interpro_accession,description > ${prefix}_intepro_summary.tsv

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
