process INTERPROSCAN {
    tag "$meta.id"
    label 'process_long'

    container "quay.io/microbiome-informatics/interproscan:5.72-103.0"

    containerOptions {
        if (workflow.containerEngine == 'singularity') {
            return "--bind ${interproscan_db}:/opt/interproscan-5.72-103.0/data"
        } else {
            return "-v ${interproscan_db}:/opt/interproscan-5.72-103.0/data"
        }
    }

    input:
    tuple val(meta), path(fasta)
    tuple path(interproscan_db), val(db_version)

    output:
    tuple val(meta), path('*.tsv') , emit: tsv
    tuple val(meta), path('*.xml') , emit: xml
    tuple val(meta), path('*.gff3'), emit: gff3
    tuple val(meta), path('*.json'), emit: json
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.name.endsWith(".gz")
    def fasta_name = fasta.name.replace(".gz", "")
    // -dp (disable precalculation) is on so no online dependency //
    """
    if ${is_compressed} ; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    interproscan.sh \\
        -cpu ${task.cpus} \\
        -i ${fasta_name} \\
        -dp \\
        ${args} \\
        --output-file-base ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        InterProScan: \$(interproscan.sh --version | grep -o "InterProScan version [0-9.-]*" | sed "s/InterProScan version //")
        InterProScan database: $db_version
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.{tsv,xml,gff3,json}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        InterProScan: \$(interproscan.sh --version | grep -o "InterProScan version [0-9.-]*" | sed "s/InterProScan version //")
        InterProScan database: $db_version
    END_VERSIONS
    """
}
