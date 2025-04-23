process MERGE_ANTISMASH_JSON {
    label 'process_medium'
    container 'ghcr.io/jqlang/jq:1.7.1'

    input:
    tuple val(meta), path(jsons)

    output:
    tuple val(meta), path('*merged.json')        , emit: merged_json
    path "versions.yml"                          , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    jq -s '
      reduce .[] as \$item (
        {
          version: null,
          input_file: "${meta.id}.fasta",
          taxon: null,
          schema: null,
          timings: {},
          records: []
        };
        {
          version: \$item.version,
          input_file: "${meta.id}.fasta",
          taxon: \$item.taxon,
          schema: \$item.schema,
          timings: (.timings + (\$item.timings // {})),
          records: (.records + (\$item.records // []))
        }
      )
    ' \$json_files > ${prefix}_merged.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jq: \$(jq -V 2>&1 | sed 's/jq-//g')
    END_VERSIONS
    """
}
