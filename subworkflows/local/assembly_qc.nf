/* LOCAL */
include { FILTER_ASSEMBLY           } from '../../modules/local/filter_assembly'
include { DECONTAMINATE_ASSEMBLIES  } from '../ebi-metagenomics/decontaminate_assemblies/main'
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

    DECONTAMINATE_ASSEMBLIES(
        FILTER_ASSEMBLY.out.fasta.ifEmpty([])
    )
    ch_versions = ch_versions.mix(DECONTAMINATE_ASSEMBLIES.out.versions)

    // Checks viability, re-compresses as bgzip, indexes, and publishes the final
    // contigs. Single ownership of _filtered_contigs.fasta.gz as a stopgap until
    // we migrate to the workflow-level outputs.
    INDEX_AND_PUBLISH_CONTIGS(
        DECONTAMINATE_ASSEMBLIES.out.cleaned_contigs
    )
    ch_versions = ch_versions.mix(INDEX_AND_PUBLISH_CONTIGS.out.versions)

    QUAST(
        ch_assembly.mix( INDEX_AND_PUBLISH_CONTIGS.out.filtered_contigs.ifEmpty([]) ).groupTuple()
    )
    ch_versions = ch_versions.mix(QUAST.out.versions)

    emit:
    assembly_qc_pass                   = INDEX_AND_PUBLISH_CONTIGS.out.filtered_contigs
    qc_failed_assemblies               = FILTER_ASSEMBLY.out.exit_reason.mix(INDEX_AND_PUBLISH_CONTIGS.out.exit_reason)
    quast_report_tsv                   = QUAST.out.tsv
    phix_contaminated_contigs_tsv      = DECONTAMINATE_ASSEMBLIES.out.phix_contaminated_contigs_tsv
    human_contaminated_contigs_tsv     = DECONTAMINATE_ASSEMBLIES.out.human_contaminated_contigs_tsv
    host_contaminated_contigs_tsv      = DECONTAMINATE_ASSEMBLIES.out.host_contaminated_contigs_tsv
    phix_contaminated_contigs_tsv_mqc  = DECONTAMINATE_ASSEMBLIES.out.phix_contaminated_contigs_tsv_mqc
    human_contaminated_contigs_tsv_mqc = DECONTAMINATE_ASSEMBLIES.out.human_contaminated_contigs_tsv_mqc
    host_contaminated_contigs_tsv_mqc  = DECONTAMINATE_ASSEMBLIES.out.host_contaminated_contigs_tsv_mqc
    versions                           = ch_versions
}
