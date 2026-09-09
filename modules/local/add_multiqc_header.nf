process ADD_MULTIQC_HEADER {
    tag "$meta.id"
    label 'process_single'

    // Reuses the gawk container purely for gzip/zcat
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
    prefix          = task.ext.prefix ?: tsv_gz.name.replaceFirst(/\.tsv\.gz$/, '')
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
