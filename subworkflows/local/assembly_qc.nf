/* LOCAL */
include { FILTER_ASSEMBLY                   } from '../../modules/local/filter_assembly'

/* NF-CORE */
include { SEQKIT_STATS as PRE_FILTER_STATS  } from '../../modules/nf-core/seqkit/stats/main'
include { SEQKIT_STATS as POST_FILTER_STATS } from '../../modules/nf-core/seqkit/stats/main'


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

    /*
    * Filter sequences based on specified criteria:
    * 1. Remove sequences shorter than the minimum contig length defined by ${params.min_contig_length}.
    * 2. Exclude sequences that contain more than 10% ambiguous bases (N).
    */
    FILTER_ASSEMBLY(
        ch_assembly
    )

    ch_versions = ch_versions.mix(FILTER_ASSEMBLY.out.versions)

    // Post filter stats //
    POST_FILTER_STATS(
        FILTER_ASSEMBLY.out.fastx
    )

    ch_versions = ch_versions.mix(FILTER_ASSEMBLY.out.versions)

    emit:
    assembly_filtered = FILTER_ASSEMBLY.out.fastx
    versions          = ch_versions
}
