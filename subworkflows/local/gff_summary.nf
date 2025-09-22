/* NF-CORE */
include { TABIX_BGZIP as TABIX_BGZIP_GFFSUMMARY } from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_GFFSUMMARY } from '../../modules/nf-core/tabix/bgziptabix/main'

/* LOCAL */
include { CREATE_GFF_SUMMARY } from '../../modules/local/create_gff_summary.nf'


workflow GFF_SUMMARY {
    take:
        ch_annotations

    main:
        ch_versions = Channel.empty()

        CREATE_GFF_SUMMARY(ch_annotations)
        ch_versions = ch_versions.mix(CREATE_GFF_SUMMARY.out.versions)

        // bgzip to get .gz and .gzi (.gzi is useful for generic paginated access to the .gff like other bgzipped .tsv files)
        TABIX_BGZIP_GFFSUMMARY(
            CREATE_GFF_SUMMARY.out.gff_summary
        )
        ch_versions = ch_versions.mix(TABIX_BGZIP_GFFSUMMARY.out.versions)

        // tabix to get .csi (.csi is like .tbi - index of genomic regions in the gff, but with better scaling for large contigs/metagenomes)
        TABIX_BGZIPTABIX_GFFSUMMARY(
            CREATE_GFF_SUMMARY.out.gff_summary  // todo: ideally tabix module would take the .gz as input instead of the uncompressed file (to avoid recompressing)
        )
        ch_versions = ch_versions.mix(TABIX_BGZIPTABIX_GFFSUMMARY.out.versions)

    emit:
        gff_summary              = CREATE_GFF_SUMMARY.out.gff_summary
        compressed_gff_summary   = TABIX_BGZIP_GFFSUMMARY.out.output
        bgzip_index              = TABIX_BGZIP_GFFSUMMARY.out.gzi
        tabix_index              = TABIX_BGZIPTABIX_GFFSUMMARY.out.gz_csi
        versions                 = ch_versions
}
