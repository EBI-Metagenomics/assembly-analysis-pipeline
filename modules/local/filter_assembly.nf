process FILTER_ASSEMBLY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/htslib_samtools_seqkit:049a7c2199a04854':
        'community.wave.seqera.io/library/htslib_samtools_seqkit:8d071a2f3053d830' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${prefix}_filtered.fasta.gz")     , emit: fasta,       optional: true
    tuple val(meta), path("${prefix}_filtered.fasta.gz.fai") , emit: fai,         optional: true
    tuple val(meta), path("${prefix}_filtered.fasta.gz.gzi") , emit: gzi,         optional: true
    tuple val(meta), env('EXIT_REASON')                      , emit: exit_reason, optional: true
    path "versions.yml"                                      , emit: versions

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
        awk '\$3 < ${params.max_contig_n_content}' > ${prefix}_nbases_filtered.tab2fx

        # Check if N-base filtering produced any sequences
        if [[ -s ${prefix}_nbases_filtered.tab2fx ]]; then
            echo "Some contigs remain after filtering..."
            echo "Converting to FASTA..."
            seqkit tab2fx ${prefix}_nbases_filtered.tab2fx \\
                --threads ${task.cpus} \\
                --out-file ${prefix}_filtered.fasta

            # bgzip to enable .gzi block index
            echo "Compressing as bgzip..."
            bgzip -@ "${task.cpus}" ${prefix}_filtered.fasta

            # using samtools as seqkit cannot index a .fasta.gz
            echo "Indexing compressed fasta..."
            samtools faidx ${prefix}_filtered.fasta.gz     # -> _filtered.fasta.gz.fai

            echo "Indexing compression archive..."
            bgzip -r ${prefix}_filtered.fasta.gz  # -> _filtered.fasta.gz.gzi
        else
            echo "No contigs after the N bases filtering"
            EXIT_REASON="insufficient_contigs_after_n_content_filtering"
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
        samtools: \$(samtools --version-only)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_filtered.fasta.gz
    touch ${prefix}_filtered.fasta.gz.fai
    touch ${prefix}_filtered.fasta.gz.gzi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
        samtools: \$(samtools --version-only)
    END_VERSIONS
    """
}
