process KEGG_ORTHOLOGS_SUMMARY {

    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/csvtk_tabix_pip_biopython:e6e033af2a05a562':
        'community.wave.seqera.io/library/csvtk_tabix_pip_biopython:7eabdb397e7420a3' }"

    input:
    tuple val(meta), path(hmmscan_concatenated_tblout)

    output:
    tuple val(meta), path("${prefix}_ko_summary.tsv.gz"),        emit: ko_summary_tsv
    tuple val(meta), path("${prefix}_ko_summary.tsv.gz.gzi"),    emit: ko_summary_tsv_gzi
    tuple val(meta), path("${prefix}_ko_per_contig.tsv.gz"),     emit: ko_per_contig_tsv
    tuple val(meta), path("${prefix}_ko_per_contig.tsv.gz.gzi"), emit: ko_per_contig_tsb_gzi
    path "versions.yml",                                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # This pipeline performs the following steps:
    # 1. Decompresses the hmmscan_concatenated_tblout file using bgzip (useful for streaming it on the website)
    # 2. Runs the hmmscan_tblout_to_tsv.py script and pipes the output
    # 3. Splits the output into two separate streams using tee:
    #    - stream 1: Extracts the first and third fields (KO ID and KO description),
    #                calculates the frequency of each unique KO ID,
    #                adds header names, reorders the fields, and saves the result to ${prefix}_ko_summary.tsv
    #    - stream 2: Extracts the first and second fields (KO ID and contig ID)
    #                and saves the result to ${prefix}_ko_per_contig.tsv - this file will be be downstream in the pipeline
    # Both TSV files are compressed with bgzip, and an index is created. The compressed files are compatible with gz, the index is
    # used on the website.

    # If there is a retry csvtk will fail with the files already exist
    rm -f ${prefix}_ko_per_contig.tsv.gz || true
    rm -f ${prefix}_ko_summary.tsv.gz || true

    hmmscan_tblout_to_tsv.py ${hmmscan_concatenated_tblout} | \\
    tee \\
        >(csvtk cut --num-cpus ${task.cpus} --tabs --no-header-row --fields 1,3 | \\
          csvtk freq --num-cpus ${task.cpus} --tabs --no-header-row --fields 1,2 --reverse --sort-by-freq | \\
          csvtk add-header --num-cpus ${task.cpus} --tabs --no-header-row --names ko,description,count | \\
          csvtk cut --num-cpus ${task.cpus} --tabs --fields count,ko,description | \\
          bgzip --stdout -@${task.cpus} --index --index-name ${prefix}_ko_summary.tsv.gz.gzi > ${prefix}_ko_summary.tsv.gz
        ) | \\
        csvtk cut --num-cpus ${task.cpus} --tabs --no-header-row --fields 1,2 | \\
        csvtk add-header --num-cpus ${task.cpus} --tabs --no-header-row --names ko,contig_id | \\
        bgzip --stdout -@${task.cpus} --index --index-name ${prefix}_ko_per_contig.tsv.gz.gzi > ${prefix}_ko_per_contig.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_ko_summary.tsv
    touch ${prefix}_ko_per_contig.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
