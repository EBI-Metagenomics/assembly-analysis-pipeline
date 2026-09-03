process FILTER_IPS_AND_FAA_BY_CONTIGS {

    tag "$meta.id"
    label 'process_medium'

    container "microbiome-informatics/seqkit_python:2.13.0_3.11"

    input:
    tuple val(meta), path(contigs_chunk), path(ips_tsv), path(faa)

    output:
    tuple val(meta), path("${prefix}_filtered_ips.tsv.gz"), path("${prefix}_filtered.faa.gz"), emit: filtered
    path "versions.yml",                                                                       emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    filter_ips_by_contigs.py \\
        --contigs ${contigs_chunk} \\
        --ips ${ips_tsv} \\
        --out ${prefix}_filtered_ips.tsv.gz | \\
    seqkit grep \\
        --threads ${task.cpus} \\
        -f - \\
        ${faa} \\
        -o ${prefix}_filtered.faa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
        python: \$( python3 --version | sed 's/Python //' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_filtered_ips.tsv.gz
    echo "" | gzip > ${prefix}_filtered.faa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
        python: \$( python3 --version | sed 's/Python //' )
    END_VERSIONS
    """
}
