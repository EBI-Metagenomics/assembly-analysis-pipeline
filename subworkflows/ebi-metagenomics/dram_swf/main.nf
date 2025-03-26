include { DRAM_DISTILL             } from '../../../modules/ebi-metagenomics/dram/distill/main'
include { SUMMARISE_FOR_DRAM_INPUT } from '../../../modules/local/summarise_for_dram_input'

workflow DRAM_SWF {

    take:
    ko_summaries       // [meta, path(tsv)]
    ko_per_contigs     // [meta, path(tsv)]
    interpro_summaries // [meta, path(tsv)]
    dbcan_overviews    // [meta, path(tsv)]

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
