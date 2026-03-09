/*
 * This process sits at the end of the end of QC (length filter → N-content
 * filter → decontamination) and serves two purposes:
 *
 * 1. Final viability check — decontamination can remove all remaining contigs
 *    even after upstream filters passed. If nothing is left the assembly is
 *    failed here with exit reason "insufficient_contigs_after_decontamination".
 *
 * 2. bgzip re-compression and indexing — This process re-compresses as bgzip and
 *    creates the .fai and .gzi indices required for random-access tools.
 *
 * It is the single publisher of the canonical _filtered_contigs.fasta.gz
 * output (which represents the most processed form of the assembly: filtered
 * and, when enabled, decontaminated). This is a stopgap until the pipeline
 * migrates to Nextflow workflow-level outputs.
 */
process INDEX_AND_PUBLISH_CONTIGS {

    tag "$meta.id"

    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/htslib_samtools_seqkit:049a7c2199a04854':
        'community.wave.seqera.io/library/htslib_samtools_seqkit:8d071a2f3053d830' }"

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("${meta.id}_filtered_contigs.fasta.gz"), emit: filtered_contigs, optional: true
    path "${meta.id}_filtered_contigs.fasta.gz.fai"              ,                         optional: true
    path "${meta.id}_filtered_contigs.fasta.gz.gzi"              ,                         optional: true
    tuple val(meta), env('EXIT_REASON')                          , emit: exit_reason,      optional: true
    path "versions.yml"                                          , emit: versions

    script:
    """
    NUM_SEQS=\$(zcat ${contigs} | grep -c '^>' || true)

    if [[ \$NUM_SEQS -eq 0 ]]; then
        EXIT_REASON="insufficient_contigs_after_decontamination"
    else
        zcat ${contigs} | bgzip -@ ${task.cpus} > ${meta.id}_filtered_contigs.fasta.gz
        samtools faidx ${meta.id}_filtered_contigs.fasta.gz
        bgzip -r ${meta.id}_filtered_contigs.fasta.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version-only)
        bgzip: \$(bgzip --version | head -1 | cut -d' ' -f3)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_filtered_contigs.fasta.gz
    touch ${meta.id}_filtered_contigs.fasta.gz.fai
    touch ${meta.id}_filtered_contigs.fasta.gz.gzi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version-only)
        bgzip: \$(bgzip --version | head -1 | cut -d' ' -f3)
    END_VERSIONS
    """
}
