process ADD_MULTIQC_HEADER {
    tag "$meta.id"
    label 'process_single'

    // Reuses the gawk container already pulled elsewhere in the decontamination
    // path (see modules/ebi-metagenomics/filterpaf) purely for gzip/zcat - no
    // awk is actually used here.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'biocontainers/gawk:5.3.0' }"

    input:
    tuple val(meta), path(tsv_gz)

    output:
    tuple val(meta), path("${prefix}_mqc.tsv.gz"), emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    // "_mqc" is appended (never reused as-is) so the output filename never matches
    // the staged input filename inside the task work dir - writing over the input
    // in-place while `zcat` is still reading it truncates the data silently. The
    // "_mqc" suffix is also what assets/multiqc_config.yml's `sp:` patterns
    // (fn: "*human_mapped_mqc.tsv" etc.) now look for, distinguishing this
    // MultiQC-only copy from the FILTERPAF output published to results
    // ("*human_mapped.tsv", no suffix - see conf/modules.config).
    prefix      = task.ext.prefix ?: tsv_gz.name.replaceFirst(/\.tsv\.gz$/, '')
    def mqc_id      = task.ext.mqc_id      ?: ''
    def mqc_section = task.ext.mqc_section ?: ''
    def mqc_title   = task.ext.mqc_title   ?: ''
    """
    {
        if [ -n "${mqc_id}" ]; then
            echo "# id: '${mqc_id}'"
            echo "# section_name: '${mqc_section}'"
            echo "# plot_type: 'table'"
            echo "# pconfig:"
            echo "#     id: '${mqc_id}_table'"
            echo "#     title: '${mqc_title}'"
        fi
        gzip -dc ${tsv_gz}
    } | gzip -c > ${prefix}_mqc.tsv.gz

    """

    stub:
    prefix = task.ext.prefix ?: tsv_gz.name.replaceFirst(/\.tsv\.gz$/, '')
    """
    touch ${prefix}_mqc.tsv.gz

    """
}
