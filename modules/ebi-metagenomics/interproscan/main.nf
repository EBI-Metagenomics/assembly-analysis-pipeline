process INTERPROSCAN {
    tag "$meta.id"
    label 'process_long'

    container 'microbiome-informatics/interproscan:5.73-104.0'

    containerOptions {
        def dataPath = "${interproscan_db}/data"
        def licensedPath = "${interproscan_db}/licensed"
        def containerArgs = []

        if (workflow.containerEngine == 'singularity') {
            containerArgs << "--bind ${dataPath}:/opt/interproscan/data"

            if (new File(licensedPath).exists()) {
                containerArgs << "--bind ${licensedPath}:/opt/interproscan/licensed"
            }

        } else {
            def workData = "${task.workDir}/${dataPath}"
            def workLicensed = "${task.workDir}/${licensedPath}"

            containerArgs << "-v ${workData}:/opt/interproscan/data"

            if (new File(workLicensed).exists()) {
                containerArgs << "-v ${workLicensed}:/opt/interproscan/licensed"
            }
        }

        return containerArgs.join(' ')
    }

    input:
    tuple val(meta), path(fasta)
    tuple path(interproscan_db), val(db_version)

    output:
    tuple val(meta), path('*.tsv.gz') , emit: tsv
    tuple val(meta), path('*.xml.gz') , emit: xml
    tuple val(meta), path('*.gff3.gz'), emit: gff3
    tuple val(meta), path('*.json.gz'), emit: json
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getExtension() == "gz"
    def fasta_file_name = fasta.name.replace(".gz", "")

    // -dp (disable precalculation) is on so no online dependency
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_file_name}
    fi

    interproscan.sh \\
        -cpu $task.cpus \\
        -i ${fasta_file_name} \\
        -dp \\
        ${args} \\
        --output-file-base ${prefix}

    gzip ${prefix}.{tsv,xml,gff3,json}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        InterProScan: \$(interproscan.sh --version | grep -o "InterProScan version [0-9.-]*" | sed "s/InterProScan version //")
        InterProScan database: $db_version
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '' > ${prefix}.{tsv,xml,gff3,json}
    gzip ${prefix}.{tsv,xml,gff3,json}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        InterProScan: \$(interproscan.sh --version | grep -o "InterProScan version [0-9.-]*" | sed "s/InterProScan version //")
        InterProScan database: $db_version
    END_VERSIONS
    """
}
