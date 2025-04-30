process GFF_SUMMARY {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    tuple val(meta), path(cds), path(ips), path(eggnog), path(dbcan_overview), path(dbcan_hmm), path(sanntis), path(antismash)

    output:
    tuple val(meta), path("*.gff.gz"), emit: gff_summary
    path "versions.yml"              , emit: versions

    script:
    """
    process_dbcan_cazys \\
        -hmm ${dbcan_hmm} \\
        -ov ${dbcan_overview} \\
        -g ${cds} \\
        -v ${params.dbcan_database_version} \\
        -o ${meta.id}_dbcan_cazys.gff

    gff_toolkit \\
        -g ${cds} \\
        -i ${ips} \\
        -e ${eggnog} \\
        -s ${sanntis} \\
        --antismash ${antismash} \\
        --dbcan-cazys ${meta.id}_dbcan_cazys.gff \\
        -o ${meta.id}_annotation_summary.gff

    gzip ${meta.id}_annotation_summary.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_annotation_summary.gff.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """
}
