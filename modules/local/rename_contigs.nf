process RENAME_CONTIGS {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyfastx:2.2.0--py39h0699b22_0' :
        'biocontainers/pyfastx:2.2.0--py39h0699b22_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*.fasta') , emit: renamed_fasta
    tuple val(meta), path('*.csv')   , emit: mapping_csv
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    rename_contigs.py --prefix ${meta.id}_ \\
    --input ${fasta} \\
    --output ${prefix}_renamed.fasta \\
    --mapping ${prefix}_mapping.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pyfastax: \$(python -c "import pyfastx; print(pyfastx.version())")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.fasta"
    touch "mapping.csv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pyfastax: \$(python -c "import pyfastx; print(pyfastx.version())")
    END_VERSIONS
    """
}
