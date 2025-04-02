process KEGG_ORTHOLOGS_SUMMARY {

    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/csvtk_tabix_pip_biopython:e6e033af2a05a562'
        : 'community.wave.seqera.io/library/csvtk_tabix_pip_biopython:7eabdb397e7420a3'}"

    input:
    tuple val(meta), path(hmmsearch_concatenated_tblout)
    path(ko_list_txt)

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
    # This mini-pipeline executes a series of steps to process the hmmsearch output:
    # 1. Decompresses the hmmsearch_concatenated_tblout file using bgzip, which is beneficial for streaming on the website.
    # 2. Extracts the contig ID from the protein ID format (e.g., transforms ERZxxx_1_2 to ERZxxx_1), including those from FGS.
    # 3. Executes the hmmsearch_tblout_to_tsv.py script to convert the hmmsearch output to TSV format and pipes the result for further processing.
    # 4. Uses the `tee` command to split the output into two distinct streams:
    #    - Stream 1:
    #      - Extracts the first (KO ID) and third (KO description) fields.
    #      - Computes the frequency of each unique KO ID.
    #      - Adds appropriate header names and reorders the fields.
    #      - Saves the processed data to a file named ${prefix}_ko_summary.tsv.
    #    - Stream 2:
    #      - Extracts the first (KO ID) and second (contig ID) fields.
    #      - Saves this data to a file named ${prefix}_ko_per_contig.tsv, which will be used in subsequent steps of the pipeline.
    # Both output TSV files are compressed using bgzip, and an index is created for each.
    # The compressed files are compatible with gz format, and the index facilitates access on the website.

    # If there is a retry csvtk will fail with the files already exist
    rm -f ${prefix}_ko_summary.tsv.gz || true
    rm -f ${prefix}_ko_summary.tsv.gz.gzi || true

    rm -f ${prefix}_ko_per_contig.tsv.gz || true
    rm -f ${prefix}_ko_per_contig.tsv.gz || true

    # Configure csvtk to use tabs
    export CSVTK_T=true

    hmmsearch_tblout_to_tsv.py ${hmmsearch_concatenated_tblout} | csvtk replace --no-header-row --fields 2 --pattern '^([A-Za-z0-9]+_[0-9]+).*\$' --replacement '\$1' | \\
    tee \\
        >(csvtk cut --num-cpus ${task.cpus} --no-header-row --fields 1 | \\
            csvtk add-header --num-cpus ${task.cpus} --no-header-row --names ko | \\
            csvtk join --num-cpus ${task.cpus} --fields "ko;knum" - ${ko_list_txt} | \\
            csvtk cut --num-cpus ${task.cpus} --fields "ko,definition" | \\
            csvtk del-header --num-cpus ${task.cpus} | \\
            csvtk freq --num-cpus ${task.cpus} --no-header-row --fields 1,2 --reverse --sort-by-freq | \\
            csvtk add-header --num-cpus ${task.cpus} --no-header-row --names ko,description,count | \\
            csvtk cut --num-cpus ${task.cpus} --fields count,ko,description | \\
            bgzip --stdout -@${task.cpus} --index --index-name ${prefix}_ko_summary.tsv.gz.gzi > ${prefix}_ko_summary.tsv.gz
        ) | \\
        csvtk cut --num-cpus ${task.cpus} --no-header-row --fields 1,2 | \\
        csvtk add-header --num-cpus ${task.cpus} --no-header-row --names ko,contig_id | \\
        bgzip --stdout -@${task.cpus} --index --index-name ${prefix}_ko_per_contig.tsv.gz.gzi > ${prefix}_ko_per_contig.tsv.gz

    # Sanity Check
    # This section verifies that the input file contains at least a few valid annotations.
    # Specifically, it checks that there are at least 2 non-comment lines in the input tblout file.
    # If the input meets this criterion, we expect the output TSV files to be non-empty.

    # Function to check if a gzipped file has at least 2 non-comment lines
    check_gz_file_lines() {
        local gz_file="\$1"

        line_count=\$(zcat "\$gz_file" | grep -v '^#' |  head -n 2 | wc -l)

        # Check if the line count is less than 2
        if [[ \$line_count -lt 2 ]]; then
            # Return false if there are less than 2 lines
            return 1
        fi
        # Return true if there are at least 2 lines
        return 0
    }

    # Check if the input file has lines and if the output files are not empty
    if check_gz_file_lines "${hmmsearch_concatenated_tblout}"; then
        # If the input file has lines, check the output files
        if ! check_gz_file_lines "${prefix}_ko_summary.tsv.gz" || ! check_gz_file_lines "${prefix}_ko_per_contig.tsv.gz"; then
            echo "Error: One or more output files do not have at least 2 lines."
            exit 1
        fi
    else
        echo "Warning: The input file does not have at least non-comments 2 lines. Skipping output file checks."
    fi

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
