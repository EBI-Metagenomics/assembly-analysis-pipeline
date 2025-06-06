/* LOCAL */
include { FILTER_ASSEMBLY } from '../../modules/local/filter_assembly'

/* NF-CORE */
include { QUAST           } from '../../modules/nf-core/quast/main'


workflow ASSEMBLY_QC {
    take:
    ch_assembly // tuple(meta, assembly_fasta)

    main:

    ch_versions = Channel.empty()

    // TODO: add the decontamination module

    /*
    * Filter sequences based on specified criteria:
    * 1. Remove sequences shorter than the minimum contig length defined by ${params.min_contig_length}.
    * 2. Exclude sequences that contain more than 10% ambiguous bases (N).
    */
    FILTER_ASSEMBLY(
        ch_assembly
    )

    ch_versions = ch_versions.mix( FILTER_ASSEMBLY.out.versions )

    QUAST(
        ch_assembly.mix( FILTER_ASSEMBLY.out.fasta.ifEmpty([]) ).groupTuple()
    )
    ch_versions = ch_versions.mix(QUAST.out.versions)

    emit:
    assembly_filtered    = FILTER_ASSEMBLY.out.fasta
    qc_failed_assemblies = FILTER_ASSEMBLY.out.exit_reason
    quast_report_tsv     = QUAST.out.tsv
    versions             = ch_versions
}
