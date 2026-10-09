process RUNDBCAN_EASYSUBSTRATE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dbcan:5.2.9--pyhdfd78af_0' :
        'quay.io/biocontainers/dbcan:5.2.9--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(input_raw_data)
    tuple val(meta2), path(input_gff), val(gff_type)
    tuple path(dbcan_db), val(dbcan_db_version)

    output:
    tuple val(meta), path("${prefix}_overview.tsv"), emit: cazyme_annotation
    tuple val(meta), path("${prefix}_dbcan_hmm_results.tsv"), emit: dbcanhmm_results
    tuple val(meta), path("${prefix}_dbcansub_hmm_results.tsv"), emit: dbcansub_results
    tuple val(meta), path("${prefix}_diamond.out"), emit: dbcandiamond_results
    tuple val(meta), path("${prefix}_cgc.gff"), emit: cgc_gff
    tuple val(meta), path("${prefix}_cgc_standard_out.tsv"), emit: cgc_standard_out
    tuple val(meta), path("${prefix}_diamond.out.tc"), emit: diamond_out_tc
    tuple val(meta), path("${prefix}_tf_hmm_results.tsv"), emit: tf_hmm_results, optional: true
    tuple val(meta), path("${prefix}_stp_hmm_results.tsv"), emit: stp_hmm_results
    tuple val(meta), path("${prefix}_total_cgc_info.tsv"), emit: total_cgc_info
    tuple val(meta), path("${prefix}_substrate_prediction.tsv"), emit: substrate_prediction
    tuple val(meta), path("${prefix}_synteny_pdf/"), optional: true, emit: synteny_pdf
    // TODO: revert this change when the migration to topics is done
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    // This are the columns of the dbcansub_hmm_results.tsv file, which will contain only a subset of the
    // headers if there are no results for the assembly
    def dbsub_output_tsv_headers = [
        'Subfam Name', 'Subfam Composition', 'Subfam EC', 'Substrate',
        'HMM Length', 'Target Name', 'Target Length', 'i-Evalue',
        'HMM From', 'HMM To', 'Target From', 'Target To', 'Coverage', 'HMM File Name'
    ].join("\t")
    """

    # Filter GFF to only include sequences present in the proteins FASTA file
    filter_gff_by_fasta_sequences.py --fasta ${input_raw_data} \\
        --gff ${input_gff} \\
        --output ${prefix}_filtered.gff

    run_dbcan easy_substrate \\
        --mode protein \\
        --db_dir ${dbcan_db} \\
        --input_raw_data ${input_raw_data} \\
        --output_dir . \\
        --input_gff ${prefix}_filtered.gff \\
        --gff_type ${gff_type} \\
        --threads ${task.cpus} \\
        ${args}

    mv overview.tsv             ${prefix}_overview.tsv
    mv dbcan_hmm_results.tsv    ${prefix}_dbcan_hmm_results.tsv
    mv dbcansub_hmm_results.tsv ${prefix}_dbcansub_hmm_results.tsv
    mv diamond.out              ${prefix}_diamond.out
    mv cgc.gff                  ${prefix}_cgc.gff
    mv cgc_standard_out.tsv     ${prefix}_cgc_standard_out.tsv
    mv diamond.out.tc           ${prefix}_diamond.out.tc
    mv stp_hmm_results.tsv      ${prefix}_stp_hmm_results.tsv
    mv total_cgc_info.tsv       ${prefix}_total_cgc_info.tsv
    mv cgc.faa                  ${prefix}_cgc.faa
    mv pul_blast.out            ${prefix}_pul_blast.out
    mv substrate_prediction.tsv ${prefix}_substrate_prediction.tsv
    mv synteny_pdf/             ${prefix}_synteny_pdf/
    if [ -f tf_hmm_results.tsv ]; then
        mv tf_hmm_results.tsv   ${prefix}_tf_hmm_results.tsv
    fi

    ##########################################################################
    # run_dbcan will produce a broken tsv if there are no results to process #
    ##########################################################################
    # The chain of warnings and errors are:
    #######################################################################
    # WARNING - No dbCAN-sub results to process
    # INFO    - Found dbcan_hmm results at results/dbCAN_hmm_results.tsv
    # ERROR   - Error loading diamond results: No columns to parse from file
    # WARNING - Missing columns in results/dbCANsub_hmm_results.tsv. Expected:
    # 'Target Name'
    # 'Subfam Name'
    # 'Subfam EC'
    # 'Target From'
    # 'Target To'
    # 'i-Evalue'
    # 
    # Found:
    # 'HMM Name'
    # 'HMM Length'
    # 'Target Name'
    # 'Target Length'
    # 'i-Evalue'
    # 'HMM From'
    # 'HMM To'
    # 'Target From'
    # 'Target To'
    # 'Coverage'
    # 'HMM File Name'
    
    # To handle this, if there is only one line in the tsv we override with all the column headers
    if [[ \$(wc -l < "${prefix}_dbcansub_hmm_results.tsv") -eq 1 ]]; then
        # I'm moving the file as otherwise I was getting a random 'cannot overwrite existing file' error
        # also, I'm not removing this file because in nfs and such systems this could be problematic as
        # there is a significant delay, so mv the file makes it easier. The cost of this is just an extra file
        # that will be deleted when the pipeline finished.
        mv ${prefix}_dbcansub_hmm_results.tsv /${prefix}_dbcansub_hmm_results.tsv_broken_headers
        echo \"${dbsub_output_tsv_headers}\" > ${prefix}_dbcansub_hmm_results.tsv
    fi
    #######################################################################

    gzip ${prefix}_*.tsv ${prefix}_*.gff ${prefix}_diamond.out.tc

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: \$(run_dbcan version | sed "s/dbCAN version: //")
        dbcan_db: "${dbcan_db_version}"
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_overview.tsv.gz
    touch ${prefix}_dbCAN_hmm_results.tsv.gz
    touch ${prefix}_dbCANsub_hmm_results.tsv.gz
    touch ${prefix}_diamond.out.gz
    touch ${prefix}_cgc.gff.gz
    touch ${prefix}_cgc_standard_out.tsv.gz
    touch ${prefix}_diamond.out.tc.gz
    touch ${prefix}_TF_hmm_results.tsv.gz
    touch ${prefix}_STP_hmm_results.tsv.gz
    touch ${prefix}_total_cgc_info.tsv.gz
    touch ${prefix}_CGC.faa.gz
    touch ${prefix}_PUL_blast.out.gz
    touch ${prefix}_substrate_prediction.tsv.gz
    mkdir -p ${prefix}_synteny_pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: \$(run_dbcan version | sed "s/dbCAN version: //")
        dbcan_db: "${dbcan_db_version}"
    END_VERSIONS
    """
}
