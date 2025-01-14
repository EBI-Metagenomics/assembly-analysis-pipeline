/* NF-CORE */
include { SEQKIT_STATS as PRE_FILTER_STATS  } from '../../modules/nf-core/seqkit/stats/main'
include { SEQKIT_STATS as POST_FILTER_STATS } from '../../modules/nf-core/seqkit/stats/main'
include { SEQKIT_SEQ                        } from '../../modules/nf-core/seqkit/seq/main'
include { SEQKIT_SPLIT2                     } from '../../modules/nf-core/seqkit/split2/main'

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

    ch_versions = ch_versions.mix(PRE_FILTER_STATS.out.versions)

    // Filter sequences shorter than ${params.min_contig_length} //
    SEQKIT_SEQ(
        ch_assembly
    )

    ch_versions = ch_versions.mix(SEQKIT_SEQ.out.versions)

    // Post filter stats //
    POST_FILTER_STATS(
        SEQKIT_SEQ.out.fastx
    )

    ch_versions = ch_versions.mix(SEQKIT_SEQ.out.versions)

    // TODO: decide when and how to chunk
    // // Split the fasta into chunks of ~max ${params.fasta_chunk_size} pb //
    // SEQKIT_SPLIT2(
    //     SEQKIT_SEQ.out.fastx
    // )
    // ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions)

    // TODO: do we need to format this?
    // Re-format the output from
    // tuple(meta, path(*.chunks))

    emit:
    assembly_filtered = SEQKIT_SEQ.out.fastx
    versions          = ch_versions
}
