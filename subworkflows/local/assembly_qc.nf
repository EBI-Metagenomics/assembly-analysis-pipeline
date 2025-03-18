/* NF-CORE */
include { SEQKIT_STATS as PRE_FILTER_STATS  } from '../../modules/nf-core/seqkit/stats/main'
include { SEQKIT_STATS as POST_FILTER_STATS } from '../../modules/nf-core/seqkit/stats/main'
include { SEQKIT_SEQ as SEQKIT_SEQ_QC       } from '../../modules/nf-core/seqkit/seq/main'

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

    // TODO: add the decontamination module

    ch_versions = ch_versions.mix(PRE_FILTER_STATS.out.versions)

    // Filter sequences shorter than ${params.min_contig_length} //
    SEQKIT_SEQ_QC(
        ch_assembly
    )

    ch_versions = ch_versions.mix(SEQKIT_SEQ_QC.out.versions)

    // Post filter stats //
    POST_FILTER_STATS(
        SEQKIT_SEQ_QC.out.fastx
    )

    ch_versions = ch_versions.mix(SEQKIT_SEQ_QC.out.versions)

    emit:
    assembly_filtered = SEQKIT_SEQ_QC.out.fastx
    versions          = ch_versions
}
