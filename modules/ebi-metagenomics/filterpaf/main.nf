process FILTERPAF {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'biocontainers/gawk:5.3.0' }"

    input:
    tuple val(meta), path(paf_file)

    output:
    tuple val(meta), path("${prefix}.txt"),           emit: mapped_contigs_txt
    tuple val(meta), path("${prefix}_mapped.tsv.gz"), emit: mapped_contigs_tsv, optional: true
    path "versions.yml",                              emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Filter PAF by query coverage and MAPQ
    awk '
        BEGIN {
            printf "sequence_id\\tquery_coverage\\tidentity\\n";
        }
        {
            query_len = \$2;
            query_start = \$3;
            query_end = \$4;
            target = \$6;
            matching_bases = \$10;
            total_bases = \$11;

            if (target == "*") {
                next;
            }

            aligned_len = query_end - query_start;
            query_coverage = aligned_len / query_len;
            identity = matching_bases / total_bases;

            if (query_coverage >= ${params.min_qcov} && identity >= ${params.min_pid}) {
                printf "%s\\t%.4f\\t%.4f\\n", \$1, query_coverage, identity;
            }
        }
    ' ${paf_file} > ${prefix}_mapped.tsv

    # Remove TSV file if it only contains the header (no contaminated contigs found)
    if [ \$(wc -l < ${prefix}_mapped.tsv) -le 1 ]; then
        rm ${prefix}_mapped.tsv
        touch ${prefix}.txt
    else
        # Extract just the sequence IDs from failed sequences
        awk '{print \$1}' ${prefix}_mapped.tsv > ${prefix}.txt
        gzip ${prefix}_mapped.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    # In stub mode, don't create the TSV file (simulating empty results)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
    END_VERSIONS
    """
}
