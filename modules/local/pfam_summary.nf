process PFAM_SUMMARY {

    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/csvtk_tabix_pip_biopython:e6e033af2a05a562':
        'community.wave.seqera.io/library/csvtk_tabix_pip_biopython:7eabdb397e7420a3' }"

    input:
    tuple val(meta), path(interproscan_tsv)

    output:
    tuple val(meta), path("${prefix}_pfam_summary.tsv.gz"),     emit: pfam_summary
    tuple val(meta), path("${prefix}_pfam_summary.tsv.gz.gzi"), emit: pfam_summary_gzi
    path "versions.yml",                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Count distinct proteins for each Pfam accession.

    # Fix quotation marks in the TSV.
    csvtk fix-quotes --tabs --no-header-row ${interproscan_tsv} |

    # Keep the Pfam rows (column 4 is "Pfam").
    csvtk filter2 --tabs --no-header-row --filter '\$4 == "Pfam"' |

    # Extract protein accession (column 1), Pfam accession (column 5)
    # and Pfam description (column 6).
    csvtk cut --tabs --no-header-row --fields 1,5,6 |

    # Keep one row per protein/Pfam accession pair.
    csvtk uniq --tabs --no-header-row --fields 1,2 |

    # Count distinct proteins per Pfam accession and description.
    csvtk freq --tabs --no-header-row --fields 2,3 --reverse --sort-by-freq |

    # Add headers and compress the TSV with bgzip, creating a .gzi index.
    csvtk add-header --tabs --no-header-row --names pfam,description,count |
    bgzip --stdout -@${task.cpus} --index \\
        --index-name ${prefix}_pfam_summary.tsv.gz.gzi \\
        > ${prefix}_pfam_summary.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_pfam_summary.tsv
    bgzip ${prefix}_pfam_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
