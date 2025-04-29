include { DRAM_DISTILL as DRAM_DISTILL_PER_ASSEMBLY    } from '../../modules/ebi-metagenomics/dram/distill/main'
include { DRAM_DISTILL as DRAM_DISTILL_PER_SAMPLESHEET } from '../../modules/ebi-metagenomics/dram/distill/main'
include { SUMMARISE_FOR_DRAM_INPUT                     } from '../../modules/local/summarise_for_dram_input'

workflow DRAM_DISTILL_SWF {

    take:
    // These should be the concatenated files - one per assembly //

    ch_proteins          // tuple [meta, path(proteins_fasta)]
    ko_per_contigs_tsv   // tuple [meta, path(tsv)]
    interproscan_tsv     // tuple [meta, path(tsv)]
    dbcan_overview       // tuple [meta, path(tsv)]

    main:

    ch_versions = Channel.empty()

    SUMMARISE_FOR_DRAM_INPUT(
        ch_proteins.join( ko_per_contigs_tsv ).join( interproscan_tsv ).join( dbcan_overview )
    )
    ch_versions = ch_versions.mix(SUMMARISE_FOR_DRAM_INPUT.out.versions)

    // Per assembly distill
    DRAM_DISTILL_PER_ASSEMBLY(
        SUMMARISE_FOR_DRAM_INPUT.out.dram_summary.groupTuple(),
        params.dram_databases
    )
    ch_versions = ch_versions.mix(DRAM_DISTILL_PER_ASSEMBLY.out.versions)

    // Whole samplesheet distill
    def samplesheet_level_dram_summary = SUMMARISE_FOR_DRAM_INPUT.out.dram_summary.map { _meta, summary -> summary }
        .collectFile(name: "samplesheet_dram_input.tsv", storeDir: "${params.outdir}/dram-distill/", keepHeader: true)

    DRAM_DISTILL_PER_SAMPLESHEET(
        samplesheet_level_dram_summary.map { summary -> [ [id:"samplesheet"], summary ] },
        params.dram_databases
    )
    ch_versions = ch_versions.mix(DRAM_DISTILL_PER_SAMPLESHEET.out.versions)

    emit:
    versions  = ch_versions
}
