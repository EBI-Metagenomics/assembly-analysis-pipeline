
/* EBI-METAGENOMICS */
include { RRNA_EXTRACTION       } from '../ebi-metagenomics/rrna_extraction/main'

workflow TAXONOMIC_ANNOTATION {

    take:
    ch_assembly // tuple (meta, assembly_fasta)

    main:

    ch_versions = Channel.empty()

    RRNA_EXTRACTION(
        ch_assembly,
        file(params.rrnas_rfam_covariance_model, checkIfExists: true),
        file(params.rrnas_rfam_claninfo, checkIfExists: true)
    )
    ch_versions = ch_versions.mix(RRNA_EXTRACTION.out.versions)

    // TODO: Sonya -- the diamond tax / CAT  bit -

    emit:
    rrna_cmsearch_deoverlap_out = RRNA_EXTRACTION.out.cmsearch_deoverlap_out
    versions                    = ch_versions
}
