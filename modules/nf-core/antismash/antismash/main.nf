process ANTISMASH_ANTISMASH {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "nf-core/antismash:8.0.0"

    input:
    tuple val(meta), path(sequence_input), path(gff)
    path databases
    val database_version

    output:
    tuple val(meta), path("${prefix}/clusterblast/*_c*.txt"), optional: true, emit: clusterblast_file
    tuple val(meta), path("${prefix}/{css,images,js}"), emit: html_accessory_files
    tuple val(meta), path("${prefix}/knownclusterblast/region*/ctg*.html"), optional: true, emit: knownclusterblast_html
    tuple val(meta), path("${prefix}/knownclusterblast/"), optional: true, emit: knownclusterblast_dir
    tuple val(meta), path("${prefix}/knownclusterblast/*_c*.txt"), optional: true, emit: knownclusterblast_txt
    tuple val(meta), path("${prefix}/svg/clusterblast*.svg"), optional: true, emit: svg_files_clusterblast
    tuple val(meta), path("${prefix}/svg/knownclusterblast*.svg"), optional: true, emit: svg_files_knownclusterblast
    tuple val(meta), path("${prefix}/${prefix}.gbk"), emit: gbk_input
    tuple val(meta), path("${prefix}/${prefix}.json"), emit: json_results
    tuple val(meta), path("${prefix}/*.log"), emit: log
    tuple val(meta), path("${prefix}/*.zip"), emit: zip
    tuple val(meta), path("${prefix}/*region*.gbk"), optional: true, emit: gbk_results
    tuple val(meta), path("${prefix}/clusterblastoutput.txt"), optional: true, emit: clusterblastoutput
    tuple val(meta), path("${prefix}/index.html"), emit: html
    tuple val(meta), path("${prefix}/knownclusterblastoutput.txt"), optional: true, emit: knownclusterblastoutput
    tuple val(meta), path("${prefix}/regions.js"), emit: json_sideloading
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.suffix ?: "${meta.id}"

    // Handle a compressed fasta file
    def is_seq_compressed = sequence_input.getExtension() == "gz" ? true : false
    def sequence_file = is_seq_compressed ? sequence_input.getBaseName() : sequence_input
    // Handle a compressed gff file
    def is_gff_compressed = gff.getExtension() == "gz" ? true : false
    def gff_file = is_gff_compressed ? gff.getBaseName() : gff

    gff_flag = gff ? "--genefinding-gff3 ${gff_file}" : ""
    """
    ## We specifically do not include on-the-fly annotations (--genefinding-tool none) as
    ## this should be run as a separate module for versioning purposes

    if [ "${is_seq_compressed}" == "true" ]; then
        # - remove in case of retry - #
        rm -f ${sequence_file}
        gzip -c -d ${sequence_input} > ${sequence_file}
    fi

    if [ "${is_gff_compressed}" == "true" ]; then
        # - remove in case of retry - #
        rm -f ${gff_file}
        gunzip -c -d ${gff} > ${gff_file}
    fi

    antismash \\
        ${args} \\
        ${gff_flag} \\
        -c ${task.cpus} \\
        --output-dir ${prefix} \\
        --output-basename ${prefix} \\
        --genefinding-tool none \\
        --logfile ${prefix}/${prefix}.log \\
        --databases ${databases} \\
        ${sequence_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash: \$(echo \$(antismash --version) | sed 's/antiSMASH //;s/-.*//g')
        antismash database: $database_version
    END_VERSIONS
    """

    stub:
    prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    """
    mkdir -p ${prefix}/css
    mkdir ${prefix}/images
    mkdir ${prefix}/js
    touch ${prefix}/NZ_CP069563.1.region001.gbk
    touch ${prefix}/NZ_CP069563.1.region002.gbk
    touch ${prefix}/css/bacteria.css
    touch ${prefix}/genome.gbk
    touch ${prefix}/genome.json
    touch ${prefix}/genome.zip
    touch ${prefix}/images/about.svg
    touch ${prefix}/index.html
    touch ${prefix}/js/antismash.js
    touch ${prefix}/js/jquery.js
    touch ${prefix}/regions.js
    touch ${prefix}/test.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash: \$(echo \$(antismash --version) | sed 's/antiSMASH //;s/-.*//g')
        antismash database: $database_version
    END_VERSIONS
    """
}
