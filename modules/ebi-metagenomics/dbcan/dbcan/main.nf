process DBCAN {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/dbcan:4.1.4--pyhdfd78af_0'
        : 'biocontainers/dbcan:4.1.4--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(gff) // only used for cluster with a proteins file
    tuple path(dbcan_db, stageAs: "dbcan_db"), val(db_version)
    val(mode)

    output:
    tuple val(meta), path("results/${prefix}_overview.txt.gz"),                             emit: overview_output
    tuple val(meta), path("results/${prefix}_*.hmm.out.gz"),                                emit: dbsub_output
    tuple val(meta), path("results/${prefix}_diamond.out.gz"),                              emit: diamond_output
    tuple val(meta), path("results/${prefix}_hmmer.out.gz"),                                emit: hmmer_output
    tuple val(meta), path("results/${prefix}_tf.out.gz"),                   optional: true, emit: tf_output
    tuple val(meta), path("results/${prefix}_tc.out.gz"),                   optional: true, emit: tc_output
    tuple val(meta), path("results/${prefix}_tp_out.gz"),                   optional: true, emit: tp_out
    tuple val(meta), path("results/${prefix}_stp.out.gz"),                  optional: true, emit: stp_out
    tuple val(meta), path("results/${prefix}_uniInput.gz"),                 optional: true, emit: uniinput
    tuple val(meta), path("results/${prefix}_cgc.out.gz"),                  optional: true, emit: cgc_output
    tuple val(meta), path("results/${prefix}_cgc.gff.gz"),                  optional: true, emit: cgc_gff
    tuple val(meta), path("results/${prefix}_cgc_standard.out.gz"),         optional: true, emit: cgc_standard_output
    tuple val(meta), path("results/${prefix}_cgc_standard.out.json.gz"),    optional: true, emit: cgc_standard_output_json
    tuple val(meta), path("results/synteny.pdf/*-syntenic.pdf.gz"),         optional: true, emit: synteny_pdfs
    tuple val(meta), path("results/${prefix}_substrate.out.gz"),            optional: true, emit: substrate_out
    path "versions.yml",                                                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    gunzip ${fasta}
    gunzip ${gff}

    run_dbcan \\
        --dia_cpu ${task.cpus} \\
        --hmm_cpu ${task.cpus} \\
        --tf_cpu ${task.cpus} \\
        --dbcan_thread ${task.cpus} \\
        --db_dir dbcan_db \\
        --out_dir results \\
        ${args} \\
        ${fasta.name.replace(".gz", "")} \\
        ${mode}

    # Bulk rename of the results, dbcan has a prefix parameter but it breaks when (for example) using --clusters and/or --cgc_substrate
    find results -type f | while read -r file; do
        mv "\$file" "\$(dirname "\$file")/${prefix}_\$(basename "\$file")"
        # rename with the prefix
        gzip "\$(dirname "\$file")/${prefix}_\$(basename "\$file")"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: 4.1.4
        dbcan database: "${db_version}"
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir results
    touch results/${prefix}_overview.txt
    touch results/${prefix}_dbout.hmm.out
    touch results/${prefix}_diamond.out
    touch results/${prefix}_hmmer.out

    gzip results/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: 4.1.4
        dbcan database: ${db_version}
    END_VERSIONS
    """
}
