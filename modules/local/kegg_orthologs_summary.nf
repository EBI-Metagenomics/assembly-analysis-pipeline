process KEGG_ORTHOLOGS_SUMMARY {

    tag "$meta.id"
    label 'process_single'

    // TODO: create image with csvtk and biopython on quay
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/csvtk_pip_biopython:251498c060995b91':
        'community.wave.seqera.io/library/csvtk_pip_biopython:e8661a0f869190e7' }"

    input:
    tuple val(meta), path(hmmscan_concatenated_tblout)

    output:
    tuple val(meta), path("${prefix}_ko_summary.tsv"),    emit: ko_summary_tsv
    tuple val(meta), path("${prefix}_ko_per_contig.tsv"), emit: ko_per_contig_tsv
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # This minipipline performs the following steps:
    # 1. Runs the hmmscan_tblout_to_tsv.py script and pipes the output
    # 2. Splits the output into two separate streams using the tee command:
    #    - stream 1: Extracts the first two fields and saves them to ko_per_contig.tsv
    #    - stream 2: Extracts the third field, calculates the frequency of each unique value in the first two fields,
    #                adds header names, reorders the fields, and saves the result to ${prefix}_ko_summary.tsv
    gunzip -c ${hmmscan_concatenated_tblout} | hmmscan_tblout_to_tsv.py | \\
        tee \
            >(csvtk cut --tabs --no-header-row --fields 1,2 > ${prefix}_ko_per_contig.tsv) \
            >(
                csvtk cut --tabs --no-header-row --fields 1,3 |
                csvtk freq --tabs --no-header-row --fields 1,2 -n |
                csvtk add-header --tabs --no-header-row --names ko,description,count |
                csvtk cut --tabs --fields count,ko,description > ${prefix}_ko_summary.tsv
            )

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_ko_summary.tsv
    touch ${prefix}_ko_per_contig.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
    END_VERSIONS
    """
}
