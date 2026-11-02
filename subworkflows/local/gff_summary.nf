/* NF-CORE */
include { HTSLIB_BGZIPTABIX as HTSLIB_BGZIPTABIX_GFFSUMMARY } from '../../modules/nf-core/htslib/bgziptabix/main'

/* LOCAL */
include { CREATE_GFF_SUMMARY } from '../../modules/local/create_gff_summary.nf'


workflow GFF_SUMMARY {
    take:
    ch_annotations

    main:
    ch_versions = channel.empty()

    CREATE_GFF_SUMMARY(ch_annotations)
    ch_versions = ch_versions.mix(CREATE_GFF_SUMMARY.out.versions)

    // bgzip to get the .gz and .gzi (.gzi is useful for generic paginated access to the .gff like other bgzipped .tsv files)
    // and tabix to get the .csi (.csi is like .tbi - index of genomic regions in the gff, but with better scaling for large contigs/metagenomes)
    HTSLIB_BGZIPTABIX_GFFSUMMARY(
        // The two empty inputs are optional: an existing .tbi/.csi index of the input (to link instead of running tabix)
        // and a regions file (tabix -R). The input is not indexed yet and we index the whole file, so neither is used.
        CREATE_GFF_SUMMARY.out.gff_summary.map { meta, gff -> [meta, gff, [], []] },
        "compress", // action
        true,       // make_index
        "gff"       // out_ext, the output is {prefix}.gff.gz
    )
    ch_versions = ch_versions.mix(HTSLIB_BGZIPTABIX_GFFSUMMARY.out.versions)

    emit:
    gff_summary              = CREATE_GFF_SUMMARY.out.gff_summary
    compressed_gff_summary   = HTSLIB_BGZIPTABIX_GFFSUMMARY.out.output
    bgzip_index              = HTSLIB_BGZIPTABIX_GFFSUMMARY.out.gzi
    tabix_index              = HTSLIB_BGZIPTABIX_GFFSUMMARY.out.index
    versions                 = ch_versions
}
