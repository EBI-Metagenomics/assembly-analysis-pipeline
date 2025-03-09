process EXTRACT_PFAM_COUNTS {

    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/csvtk:0.31.0--h9ee0642_0':
        'biocontainers/csvtk:0.31.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(interproscan_tsv)

    output:
    tuple val(meta), path("${prefix}_pfam.tsv"), emit: output
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # This csvtk minipipeline extracts Pfam rows:
    # 1. Picks the Pfam rows.
    # 2. Extracts the 5th (Pfam accession) and 6th (description) columns.
    # 3. Counts the frequency of the pfams (we use pfam accession and description because we need the description too)
    # 4. Adds headers ('pfam', 'description' and 'count') the TSV.
    # 5. Inverts the columns - we need count to the be first

    csvtk filter2 --tabs --no-header-row --filter '\$4 == "Pfam"' ${interproscan_tsv} | \\
    csvtk cut --tabs --no-header-row --fields 5,6 | \\
    csvtk freq --tabs --no-header-row --fields 1,2 -n | \\
    csvtk add-header --tabs --no-header-row --names pfam,description,count | \\
    csvtk cut --tabs --fields count,pfam,description > ${prefix}_pfam.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_pfam.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
    END_VERSIONS
    """
}
