process ANTISMASH_ANTISMASHLITE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/microbiome-informatics/antismash:7.1.0.1_2'

    input:
    tuple val(meta), path(sequence_input)
    path(databases)
    path(gff)

    output:
    tuple val(meta), path("${prefix}/clusterblast/*_c*.txt")                 , optional: true, emit: clusterblast_file
    tuple val(meta), path("${prefix}/{css,images,js}")                       , emit: html_accessory_files
    tuple val(meta), path("${prefix}/knownclusterblast/region*/ctg*.html")   , optional: true, emit: knownclusterblast_html
    tuple val(meta), path("${prefix}/knownclusterblast/")                    , optional: true, emit: knownclusterblast_dir
    tuple val(meta), path("${prefix}/knownclusterblast/*_c*.txt")            , optional: true, emit: knownclusterblast_txt
    tuple val(meta), path("${prefix}/svg/clusterblast*.svg")                 , optional: true, emit: svg_files_clusterblast
    tuple val(meta), path("${prefix}/svg/knownclusterblast*.svg")            , optional: true, emit: svg_files_knownclusterblast
    tuple val(meta), path("${prefix}/*.gbk")                                 , emit: gbk_input
    tuple val(meta), path("${prefix}/*.json")                                , emit: json_results
    tuple val(meta), path("${prefix}/*.log")                                 , emit: log
    tuple val(meta), path("${prefix}/*.zip")                                 , emit: zip
    tuple val(meta), path("${prefix}/*region*.gbk")                          , optional: true, emit: gbk_results
    tuple val(meta), path("${prefix}/clusterblastoutput.txt")                , optional: true, emit: clusterblastoutput
    tuple val(meta), path("${prefix}/index.html")                            , emit: html
    tuple val(meta), path("${prefix}/knownclusterblastoutput.txt")           , optional: true, emit: knownclusterblastoutput
    tuple val(meta), path("${prefix}/regions.js")                            , emit: json_sideloading
    path "versions.yml"                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def gff_flag = gff ? "--genefinding-gff3 ${gff}" : ""
    def is_compressed = sequence_input.getExtension() == "gz" ? true : false
    def sequence_file = is_compressed ? sequence_input.getBaseName() : sequence_input
    """
    ## We specifically do not include on-the-fly annotations (--genefinding-tool none) as
    ## this should be run as a separate module for versioning purposes    
    if [ "${is_compressed}" == "true" ]; then
        # - remove in case of retry - #
        rm -f ${sequence_file}
        gzip -c -d ${sequence_input} > ${sequence_file}
    fi

    #######################
    # Clean the GFF3 file #
    #######################
    # TODO: This is a patch to allow antiSMASH to ingest the InterProScan GFF3. It should probably be a Python script.
    # InterProScan runs on the proteins, which are labeled based on the Combined Gene Caller conventions (which come from Prodigal and FragGeneScan).
    # AntiSMASH needs the ID (the first column) to match the contig, which is the first element of this AWK expression, to convert the following:
    # A FGS example: contig-2-length-2192_210_335_+ -> contig-2-length-2192 (the contig name).
    # The second part is to change the third column from IPS to CDS (which is what AS expects... I know).
    # Finally, the last part is to remove the Name="<name>" from IPS, as AS doesn't like repeated values there... again, I know.
    # Example: Name=mobidb-lite will be in every prediction that comes from MobiDBLite.

    # - remove in case of retry - #
    rm -f input.gff
    cat ${gff} | awk 'BEGIN{FS=OFS="\t"} /^#/{print; next} {split(\$1,a,"_"); \$1=a[1]; \$3="CDS"; sub(/Name=[^;]+;/, "", \$9); print}' > input.gff

    antismash \\
        ${args} \\
        --genefinding-gff3 input.gff \\
        -c ${task.cpus} \\
        --output-dir ${prefix} \\
        --output-basename ${prefix} \\
        --genefinding-tool none \\
        --logfile $prefix/${prefix}.log \\
        --databases $databases \\
        $sequence_input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash-lite: \$(echo \$(antismash --version) | sed 's/antiSMASH //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def VERSION = '7.1.0' // WARN: Version information not provided by tool during stub run. Please update this string when bumping container versions.
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
        antismash-lite: $VERSION
    END_VERSIONS
    """
}
