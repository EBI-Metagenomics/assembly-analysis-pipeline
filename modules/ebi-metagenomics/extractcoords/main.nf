
process EXTRACTCOORDS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    tuple val(meta), path(easel_coords_fasta)
    tuple val(meta2), path(matched_seqs_with_coords)

    output:
    tuple val(meta), path("sequence-categorisation/*SSU.fasta.gz")        , optional: true, emit: ssu_fasta
    tuple val(meta), path("sequence-categorisation/*LSU.fasta.gz")        , optional: true, emit: lsu_fasta
    tuple val(meta), path("sequence-categorisation/*5S.fasta.gz")         , optional: true, emit: fiveS_fasta
    tuple val(meta), path("sequence-categorisation/*5_8S.fasta.gz")       , optional: true, emit: five_eightS_fasta
    tuple val(meta), path("sequence-categorisation/*other_ncRNA.fasta.gz"), optional: true, emit: ncrna_fasta
    tuple val(meta), path("*concat_SSU_LSU_coords.txt")                   , emit: concat_ssu_lsu_coords
    path "versions.yml"                                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    get_subunits -i $easel_coords_fasta -n ${prefix}

    get_subunits_coords -i $matched_seqs_with_coords -s SSU -l LSU

    cat SSU_coords LSU_coords > ${prefix}_concat_SSU_LSU_coords.txt

    gzip sequence-categorisation/*.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_concat_SSU_LSU_coords.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """
}
