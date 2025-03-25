include { DRAM                  } from '../../../modules/ebi-metagenomics/dram/distill/main'
include { SUMMARISEFORDRAMINPUT } from '../../../modules/local/summarisefordraminput/main'

workflow DRAM_SWF {

    take:
    root_path     //

    main:

    ch_versions = Channel.empty()

    SUMMARISEFORDRAMINPUT (
        root_path
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
