/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC as MULTIQC_PER_ASSEMBLY    } from '../modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_PER_SAMPLESHEET } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                   } from 'plugin/nf-schema'
include { paramsSummaryMultiqc               } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML             } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText             } from '../subworkflows/local/utils_nfcore_assembly_analysis_pipeline_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   NF-CORE MODULES and SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SEQKIT_SPLIT2                     } from '../modules/nf-core/seqkit/split2/main'
include { PIGZ_UNCOMPRESS as PIGZ_CONTIGS   } from '../modules/nf-core/pigz/uncompress/main'
include { PIGZ_UNCOMPRESS as PIGZ_PROTEINS  } from '../modules/nf-core/pigz/uncompress/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   EBI-METAGENOMICS MODULES and SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { ASSEMBLY_QC                       } from '../subworkflows/local/assembly_qc'
include { COMBINED_GENE_CALLER              } from '../subworkflows/ebi-metagenomics/combined_gene_caller/main'
include { CONTIGS_TAXONOMIC_CLASSIFICATION  } from '../subworkflows/ebi-metagenomics/contigs_taxonomic_classification/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   LOCAL MODULES and SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RENAME_CONTIGS         } from '../modules/local/rename_contigs'
include { RNA_ANNOTATION         } from '../subworkflows/local/rna_annotation'
include { FUNCTIONAL_ANNOTATION  } from '../subworkflows/local/functional_annotation'
include { PATHWAYS_AND_SYSTEMS   } from '../subworkflows/local/pathways_and_systems'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ASSEMBLY_ANALYSIS_PIPELINE {
    take:
    ch_assembly // tuple( meta, assembly_fasta )

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    /*
     * Rename the contigs using the prodivded prefix, seqs will be named >{prefix}_{n}
     * there n is just an autoincrement
    */
    RENAME_CONTIGS(
        ch_assembly
    )

    /*
    * The first step is to:
    * - Gather some statistics about the assembly
    * - Filter by length
    * - TODO: Remove human and any other host specified - https://www.ebi.ac.uk/panda/jira/browse/EMG-7487
    */
    ASSEMBLY_QC(
        RENAME_CONTIGS.out.renamed_fasta
    )
    ch_versions = ch_versions.mix(ASSEMBLY_QC.out.versions)

    /*
     * Run the RNA detection subworklow
    */
    RNA_ANNOTATION(
        ASSEMBLY_QC.out.assembly_filtered
    )

    /*
    * Protein prediction with the combined-gene-caller, and masking the rRNAs genes
    */
    // We need to sync the sequences and the rRNA outputs //
    ASSEMBLY_QC.out.assembly_filtered.join(RNA_ANNOTATION.out.ssu_lsu_coords).multiMap { meta, assembly_fasta, ssu_lsu_coords ->
        assembly: [meta, assembly_fasta]
        ssu_lsu_coords: [meta, ssu_lsu_coords]
    }.set {
        ch_cgc
    }

    // TODO: handle LR - FGS flip parameter //
    COMBINED_GENE_CALLER(
        ch_cgc.assembly,
        ch_cgc.ssu_lsu_coords // used to mask the rRNA genes in the assembly
    )
    ch_versions = ch_versions.mix(COMBINED_GENE_CALLER.out.versions)

    /*
     * Taxonomic classification of the contigs with CATPACK
     * This outside of the taxonomical_annotation suboworkflow because it has a dependency with the
     * CGC predicted proteins
    */
    // TOOD: handle gzip files in this subworkflow instead of uncompress this files
    CONTIGS_TAXONOMIC_CLASSIFICATION(
        PIGZ_CONTIGS(ASSEMBLY_QC.out.assembly_filtered).file,
        PIGZ_PROTEINS(COMBINED_GENE_CALLER.out.faa).file,
        [[id: "cat_diamond_db"], file(params.cat_diamond_database, checkIfExists: true)],
        [[id: "cat_taxonomy_db"], file(params.cat_taxonomy_database, checkIfExists: true)]
    )
    ch_versions = ch_versions.mix(CONTIGS_TAXONOMIC_CLASSIFICATION.out.versions)

    /*********************/
    /* PROTEINS CHUNKING */
    /*********************/

    // Chunk the fasta into files with at most >= params
    SEQKIT_SPLIT2(
        COMBINED_GENE_CALLER.out.faa,
        params.proteins_chunksize,
    )
    ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions)

    def ch_protein_chunks = SEQKIT_SPLIT2.out.assembly.transpose()

    /*
    * Annotation of the proteins.
    */
    FUNCTIONAL_ANNOTATION(
        COMBINED_GENE_CALLER.out.faa.join( COMBINED_GENE_CALLER.out.gff ),
        ch_protein_chunks
    )
    ch_versions = ch_versions.mix(FUNCTIONAL_ANNOTATION.out.versions)

    // TODO: print each of the following ones to see how to handle meta, study, analyses
    //def ko_summaries = FUNCTIONAL_ANNOTATION.out.ko_summary_tsv.collect()
    //def ko_per_contigs = FUNCTIONAL_ANNOTATION.out.ko_per_contig_tsv.collect()
    //def interpro_summaries = FUNCTIONAL_ANNOTATION.out.interproscan_tsv.collect()
    //def dbcan_overviews = FUNCTIONAL_ANNOTATION.out.dbcan_overview.collect()

    //DRAM_SWF(
    //    FUNCTIONAL_ANNOTATION.out.ko_summary_tsv.collect(),
    //    FUNCTIONAL_ANNOTATION.out.ko_per_contig_tsv.collect(),
    //    FUNCTIONAL_ANNOTATION.out.interproscan_tsv.collect(),
    //    FUNCTIONAL_ANNOTATION.out.dbcan_overview.collect()
    //)

    /*
    * Pathway and systems annotations
    */
    PATHWAYS_AND_SYSTEMS(
        ASSEMBLY_QC.out.assembly_filtered.join(
            COMBINED_GENE_CALLER.out.faa
        ).join(
            COMBINED_GENE_CALLER.out.gff
        ).join(
            FUNCTIONAL_ANNOTATION.out.interproscan_tsv
        ),
        FUNCTIONAL_ANNOTATION.out.kegg_orthologs_summary_tsv,
        ch_protein_chunks
    )
    ch_versions = ch_versions.mix(PATHWAYS_AND_SYSTEMS.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'assembly_analysis_pipeline_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = Channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : Channel.fromPath("${projectDir}/assets/mgnify_wordmark_dark_on_light.png", checkIfExists: true)

    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )
    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    /**************************************/
    /* MultiQC report for the samplesheet */
    /**************************************/

    common_files = ch_multiqc_files.collect()

    MULTIQC_PER_ASSEMBLY(
        ASSEMBLY_QC.out.quast_report_tsv,
        common_files,
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    MULTIQC_PER_SAMPLESHEET(
        ASSEMBLY_QC.out.quast_report_tsv.collect(),
        common_files,
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    multiqc_report = MULTIQC_PER_SAMPLESHEET.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                                // channel: [ path(versions.yml) ]
}
