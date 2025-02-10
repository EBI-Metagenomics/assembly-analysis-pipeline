/* NF-CORE */
include { SEQKIT_STATS as PRE_FILTER_STATS  } from '../../modules/nf-core/seqkit/stats/main'
include { SEQKIT_STATS as POST_FILTER_STATS } from '../../modules/nf-core/seqkit/stats/main'
include { SEQKIT_SEQ                        } from '../../modules/nf-core/seqkit/seq/main'

/* EBI-METAGENOMICS */


workflow ASSEMBLY_QC {
    take:
    ch_assembly // tuple(meta, assembly_fasta)

    main:

    ch_versions = Channel.empty()

    // Pre filtering stats //
    PRE_FILTER_STATS(
        ch_assembly
    )


    // TODOs:
    // Check 5 trimming, ambigous bases and other bits
    // Decont for human stuff

    ch_versions = ch_versions.mix(PRE_FILTER_STATS.out.versions)

    // Filter sequences shorter than ${params.min_contig_length} //
    // TODO: metaT size is 200 for MAP - 1kb
    SEQKIT_SEQ(
        ch_assembly
    )

    ch_versions = ch_versions.mix(SEQKIT_SEQ.out.versions)

    // Post filter stats //
    POST_FILTER_STATS(
        SEQKIT_SEQ.out.fastx
    )

    ch_versions = ch_versions.mix(SEQKIT_SEQ.out.versions)

    emit:
    assembly_filtered = SEQKIT_SEQ.out.fastx
    versions          = ch_versions
}
