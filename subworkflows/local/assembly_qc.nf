/* LOCAL */
include { FILTER_ASSEMBLY           } from '../../modules/local/filter_assembly'
include { ASSEMBLY_DECONTAMINATION  } from '../ebi-metagenomics/assembly_decontamination/main'
include { INDEX_AND_PUBLISH_CONTIGS } from '../../modules/local/index_and_publish_contigs'

/* NF-CORE */
include { QUAST                    } from '../../modules/nf-core/quast/main'


workflow ASSEMBLY_QC {
    take:
    ch_assembly // tuple(meta, assembly_fasta)

    main:

    ch_versions = channel.empty()

    /*
    * Filter sequences based on specified criteria:
    * 1. Remove sequences shorter than the minimum contig length defined by ${params.min_contig_length}.
    * 2. Exclude sequences that contain more than 10% ambiguous bases (N).
    * 3. Run the contigs decontamination module (a. and b. are optional)
    *    a. To remove human and phyx containamted contigs
    *.   b. Contigs contaninated with whichever host genome was selected
    */

    FILTER_ASSEMBLY(
        ch_assembly
    )

    ch_versions = ch_versions.mix( FILTER_ASSEMBLY.out.versions )

    ASSEMBLY_DECONTAMINATION(
        FILTER_ASSEMBLY.out.fasta.ifEmpty([])
    )
    ch_versions = ch_versions.mix(ASSEMBLY_DECONTAMINATION.out.versions)

    // Checks viability, re-compresses as bgzip, indexes, and publishes the final
    // contigs. Single ownership of _filtered_contigs.fasta.gz as a stopgap until
    // we migrate to the workflow-level outputs.
    INDEX_AND_PUBLISH_CONTIGS(
        ASSEMBLY_DECONTAMINATION.out.cleaned_contigs
    )
    ch_versions = ch_versions.mix(INDEX_AND_PUBLISH_CONTIGS.out.versions)

    QUAST(
        ch_assembly.mix( INDEX_AND_PUBLISH_CONTIGS.out.filtered_contigs.ifEmpty([]) ).groupTuple()
    )
    ch_versions = ch_versions.mix(QUAST.out.versions)

    emit:
    assembly_qc_pass                = INDEX_AND_PUBLISH_CONTIGS.out.filtered_contigs
    qc_failed_assemblies            = FILTER_ASSEMBLY.out.exit_reason.mix(INDEX_AND_PUBLISH_CONTIGS.out.exit_reason)
    quast_report_tsv                = QUAST.out.tsv
    phix_contaminated_contigs_tsv   = ASSEMBLY_DECONTAMINATION.out.phix_contaminated_contigs_tsv
    human_contaminated_contigs_tsv  = ASSEMBLY_DECONTAMINATION.out.human_contaminated_contigs_tsv
    host_contaminated_contigs_tsv   = ASSEMBLY_DECONTAMINATION.out.host_contaminated_contigs_tsv
    versions                        = ch_versions
}
