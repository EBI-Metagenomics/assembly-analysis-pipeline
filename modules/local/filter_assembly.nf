process FILTER_ASSEMBLY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.8.1--h9ee0642_0':
        'biocontainers/seqkit:2.8.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${prefix}_filtered.fasta.gz") , emit: fasta,       optional: true
    tuple val(meta), env('EXIT_REASON')                  , emit: exit_reason, optional: true
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix  = task.ext.prefix ?: "${meta.id}"
    """
    seqkit \\
        seq \\
        --min-len ${params.min_contig_length} \\
        --threads $task.cpus \\
        $assembly > ${prefix}_len_filtered.fasta

    # Check if length filtering produced any sequences
    if [[ \$(wc -l < ${prefix}_len_filtered.fasta) -eq 0 ]]; then
        echo "No contigs after the length filtering"
        EXIT_REASON="insufficient_contigs_after_length_filtering"
    else
        echo "Filtering sequences by N-base content (< 10%)..."

        seqkit fx2tab ${prefix}_len_filtered.fasta \\
            --threads ${task.cpus} \\
            --base-content N | \\
        awk '\$3 < 10' > ${prefix}_nbases_filtered.tab2fx

        # Check if N-base filtering produced any sequences
        if [[ -s ${prefix}_nbases_filtered.tab2fx ]]; then
            seqkit tab2fx ${prefix}_nbases_filtered.tab2fx \\
                --threads ${task.cpus} \\
                --out-file ${prefix}_filtered.fasta.gz
        else
            echo "No contigs after the N bases filtering"
            EXIT_REASON="insufficient_contigs_after_n_bases_filtering"
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_filtered.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """
}
