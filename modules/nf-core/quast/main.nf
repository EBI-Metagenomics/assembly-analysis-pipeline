process QUAST {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quast:5.2.0--py39pl5321heaaa4ec_4' :
        'biocontainers/quast:5.2.0--py39pl5321heaaa4ec_4' }"

    input:
    tuple val(meta) , path(consensus)

    output:
    tuple val(meta), path("${prefix}")                   , emit: results
    tuple val(meta), path("${prefix}.tsv.gz")            , emit: tsv
    tuple val(meta), path("${prefix}_transcriptome.tsv") , optional: true , emit: transcriptome
    tuple val(meta), path("${prefix}_misassemblies.tsv") , optional: true , emit: misassemblies
    tuple val(meta), path("${prefix}_unaligned.tsv")     , optional: true , emit: unaligned
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args   ?: ''
    prefix        = task.ext.prefix ?: "${meta.id}"
    """
    metaquast.py \\
        --output-dir $prefix \\
        --threads $task.cpus \\
        $args \\
        ${consensus.join(' ')}
    
    gzip ${prefix}/report.tsv

    ln -s ${prefix}/report.tsv.gz ${prefix}.tsv.gz
    [ -f  ${prefix}/contigs_reports/all_alignments_transcriptome.tsv ] && ln -s ${prefix}/contigs_reports/all_alignments_transcriptome.tsv ${prefix}_transcriptome.tsv
    [ -f  ${prefix}/contigs_reports/misassemblies_report.tsv         ] && ln -s ${prefix}/contigs_reports/misassemblies_report.tsv ${prefix}_misassemblies.tsv
    [ -f  ${prefix}/contigs_reports/unaligned_report.tsv             ] && ln -s ${prefix}/contigs_reports/unaligned_report.tsv ${prefix}_unaligned.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix        = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p $prefix
    touch $prefix/report.tsv
    gzip ${prefix}.tsv
    touch $prefix/report.html
    touch $prefix/report.pdf
    touch $prefix/quast.log
    touch $prefix/transposed_report.txt
    touch $prefix/transposed_report.tex
    touch $prefix/icarus.html
    touch $prefix/report.tex
    touch $prefix/report.txt

    mkdir -p $prefix/basic_stats
    touch $prefix/basic_stats/cumulative_plot.pdf
    touch $prefix/basic_stats/Nx_plot.pdf
    touch $prefix/basic_stats/genome_GC_content_plot.pdf
    touch $prefix/basic_stats/GC_content_plot.pdf

    mkdir -p $prefix/icarus_viewers
    touch $prefix/icarus_viewers/contig_size_viewer.html

    ln -s $prefix/report.tsv ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """
}
