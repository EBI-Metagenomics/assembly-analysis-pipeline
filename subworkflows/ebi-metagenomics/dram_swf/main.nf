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

    SUMMARISE_FOR_DRAM_INPUT (
        ko_summaries,
        ko_per_contigs,
        interpro_summaries,
        dbcan_overviews
    )
    ch_versions = ch_versions.mix(SUMMARISE_FOR_DRAM_INPUT.out.versions)

    DRAM_DISTILL(
        SUMMARISE_FOR_DRAM_INPUT.out.dram_summary, 
        params.dram_databases
    )
    ch_versions = ch_versions.mix(DRAM_DISTILL.out.versions)

    emit:
    dram_tsv  = DRAM_DISTILL.out.out_tsv
    dram_html = DRAM_DISTILL.out.html
    versions  = ch_versions
}
