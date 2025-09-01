/* NF-CORE */
include { TABIX_BGZIP as TABIX_BGZIP_GFFSUMMARY            } from '../../modules/nf-core/tabix/bgzip/main'

/* LOCAL */
include { CREATE_GFF_SUMMARY } from '../../modules/local/create_gff_summary.nf'


workflow GFF_SUMMARY {
    take:
        ch_annotations // tuple (meta, path(cds), path(ips), path(eggnog), path(dbcan_overview), path(dbcan_hmm), path(sanntis), path(antismash))

    main:
        ch_versions = Channel.empty()

        CREATE_GFF_SUMMARY(ch_annotations)
        ch_versions = ch_versions.mix(CREATE_GFF_SUMMARY.out.versions)

        TABIX_BGZIP_GFFSUMMARY(
            CREATE_GFF_SUMMARY.out.gff_summary
        )
        ch_versions = ch_versions.mix(TABIX_BGZIP_GFFSUMMARY.out.versions)

    emit:
        gff_summary  = CREATE_GFF_SUMMARY.out.gff_summary
        compressed_gff_summary  = TABIX_BGZIP_GFFSUMMARY.out.output
        tabix_index   = TABIX_BGZIP_GFFSUMMARY.out.gzi
        versions    = ch_versions
}
