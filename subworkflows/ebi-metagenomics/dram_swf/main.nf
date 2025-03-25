include { DRAM                  } from '../../../modules/ebi-metagenomics/dram/distill/main'
include { SUMMARISEFORDRAMINPUT } from '../../../modules/local/summarisefordraminput/main'

workflow DRAM_SWF {

    take:
    tuple val(meta), path(ko_summaries)
    tuple val(meta), path(ko_per_contigs)
    tuple val(meta), path(interpro_summaries)
    tuple val(meta), path(dbcan_overviews)

    main:

    ch_versions = Channel.empty()

    SUMMARISEFORDRAMINPUT (
        ko_summaries,
        ko_per_contigs,
        interpro_summaries,
        dbcan_overviews
    )
    ch_versions = ch_versions.mix(SUMMARISEFORDRAMINPUT.out.versions)

    DRAM(
        SUMMARISEFORDRAMINPUT.out.tsv
    )
    ch_versions = ch_versions.mix(DRAM.out.versions)

    emit:
    dram_tsv  = DRAM.out.tsv	
    dram_html = DRAM.out.html
    versions  = ch_versions
}
