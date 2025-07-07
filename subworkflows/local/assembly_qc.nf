/* LOCAL */
include { FILTER_ASSEMBLY          } from '../../modules/local/filter_assembly'
include { ASSEMBLY_DECONTAMINATION } from '../ebi-metagenomics/assembly_decontamination/main'

/* NF-CORE */
include { QUAST                    } from '../../modules/nf-core/quast/main'


workflow ASSEMBLY_QC {
    take:
    ch_assembly // tuple(meta, assembly_fasta)

    main:

    ch_versions = Channel.empty()

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

    QUAST(
        ch_assembly.mix( ASSEMBLY_DECONTAMINATION.out.cleaned_contigs.ifEmpty([]) ).groupTuple()
    )
    ch_versions = ch_versions.mix(QUAST.out.versions)

    emit:
    assembly_qc_pass     = ASSEMBLY_DECONTAMINATION.out.cleaned_contigs
    qc_failed_assemblies = FILTER_ASSEMBLY.out.exit_reason
    quast_report_tsv     = QUAST.out.tsv
    versions             = ch_versions
}
