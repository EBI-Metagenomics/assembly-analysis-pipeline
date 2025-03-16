process DIAMOND_RHEACHEBI {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/diamond_mgnify-pipelines-toolkit:059c7a2806a8307e' :
        'community.wave.seqera.io/library/diamond_mgnify-pipelines-toolkit:d550c08ee965686f' }"

    input:
    tuple val(meta) , path(fasta)
    path(db)
    path(rhea2chebi)
    val blast_columns

    output:
    tuple val(meta), path("*_rhea2proteins.tsv.gz"), emit: rhea2proteins_tsv
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def columns = blast_columns ? "${blast_columns}" : ''
    def outfmt = '6' // txt
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    diamond \\
        blastp \\
        --threads ${task.cpus} \\
        --evalue 1e-10 \\
        --query-cover 80 \\
        --subject-cover 80 \\
        --id 50 \\
        -k 0 \\
        --db ${db} \\
        --query ${fasta} \\
        --outfmt ${outfmt} ${columns} | \\
    add_rhea_chebi_annotation_patched.py \\
        --diamond_hits - \\
        --proteins ${fasta_name} \\
        --rhea2chebi ${rhea2chebi} \\
        --output ${prefix}_rhea2proteins.tsv

    gzip ${prefix}_rhea2proteins.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def out_ext = "txt"
    """
    touch ${prefix}.${out_ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """
}
