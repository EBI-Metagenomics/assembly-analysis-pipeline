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
include { SEQKIT_SPLIT2                    } from '../modules/nf-core/seqkit/split2/main'
include { PIGZ_UNCOMPRESS as PIGZ_CONTIGS  } from '../modules/nf-core/pigz/uncompress/main'
include { PIGZ_UNCOMPRESS as PIGZ_PROTEINS } from '../modules/nf-core/pigz/uncompress/main'
include { FIND_UNPIGZ                      } from '../modules/nf-core/find/unpigz/main'
include { GT_GFF3VALIDATOR                 } from '../modules/nf-core/gt/gff3validator/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EBI-METAGENOMICS MODULES and SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { DOWNLOAD_FROM_FIRE                 } from '../modules/ebi-metagenomics/download_from_fire/main'
include { ASSEMBLY_QC                        } from '../subworkflows/local/assembly_qc'
include { COMBINED_GENE_CALLER               } from '../subworkflows/ebi-metagenomics/combined_gene_caller/main'
include { CONTIGS_TAXONOMIC_CLASSIFICATION   } from '../subworkflows/ebi-metagenomics/contigs_taxonomic_classification/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL MODULES and SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RENAME_CONTIGS                     } from '../modules/local/rename_contigs'
include { GFF_SUMMARY                        } from '../subworkflows/local/gff_summary'
include { RNA_ANNOTATION                     } from '../subworkflows/local/rna_annotation'
include { FUNCTIONAL_ANNOTATION              } from '../subworkflows/local/functional_annotation'
include { PATHWAYS_AND_SYSTEMS               } from '../subworkflows/local/pathways_and_systems'

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
     * FIRE download: Download files from EBI FIRE system if enabled
     * This is only available on EBI network
     */
    if ( params.use_fire_download ) {
        // This module will transform the FTP links to FIRE-S3 ones and it will download them
        DOWNLOAD_FROM_FIRE(
            ch_assembly.map { meta, assembly_fasta -> {
                    [meta, [assembly_fasta]] // This is workaround as the module expects a list of accessions (for fwd/rev reads)
                }
            }
        )
        ch_versions = ch_versions.mix(DOWNLOAD_FROM_FIRE.out.versions)

        // Replace original assembly paths with downloaded files
        ch_assembly = ch_assembly
            .join(DOWNLOAD_FROM_FIRE.out.downloaded_files)
            .map { meta, _assembly_fasta, downloaded_files ->
                [meta, downloaded_files] // It's just the one file
            }
    }

    /*
     * Rename the contigs using the provided prefix, seqs will be named >{prefix}_{n}
     * there n is just an autoincrement
    */
    RENAME_CONTIGS(
        ch_assembly
    )
    ch_versions = ch_versions.mix(RENAME_CONTIGS.out.versions)

    /*
    * The first step is to:
    * - Gather some statistics about the assembly
    * - Filter by length
    * - Remove human, phiX and/or host/contanimation data from the contigs.
    */

    /*
    We need to adjust the meta, the human and PhiX references can be either set in the samplesheet, otherwise
    they will be taken from the parameters
    */
    def renamed_contigs = RENAME_CONTIGS.out.renamed_fasta.map {
        meta, contigs -> {
            def new_meta = meta.clone()
            if ( params.skip_decontamination ) {
                // We remove the refs from the meta to skip that code
                new_meta.contaminant_reference = null
                new_meta.human_reference = null
                new_meta.phix_reference  = null
            } else {
                // If the user didn't provide refs, use the ones in the parameters
                if (!meta.contaminant_reference) {
                    new_meta.contaminant_reference = params.contaminant_reference ?: null
                }
                if (!meta.human_reference) {
                    new_meta.human_reference = params.human_reference ?: null
                }
                if (!meta.phix_reference) {
                    new_meta.phix_reference = params.phix_reference ?: null
                }
            }
            [new_meta, contigs]
        }
    }

    ASSEMBLY_QC(
        renamed_contigs
    )
    ch_versions = ch_versions.mix(ASSEMBLY_QC.out.versions)

    /*
     * Run the RNA detection subworklow
    */
    RNA_ANNOTATION(
        ASSEMBLY_QC.out.assembly_qc_pass
    )
    ch_versions = ch_versions.mix(RNA_ANNOTATION.out.versions)

    /*
    * Protein prediction with the combined-gene-caller, and masking the rRNAs genes
    */
    // We need to sync the sequences and the rRNA outputs //
    ASSEMBLY_QC.out.assembly_qc_pass
        .join(RNA_ANNOTATION.out.ssu_lsu_coords)
        .multiMap { meta, assembly_fasta, ssu_lsu_coords ->
            assembly: [meta, assembly_fasta]
            ssu_lsu_coords: [meta, ssu_lsu_coords]
        }
        .set {
            ch_cgc
        }

    COMBINED_GENE_CALLER(
        ch_cgc.assembly,
        ch_cgc.ssu_lsu_coords,
    )
    ch_versions = ch_versions.mix(COMBINED_GENE_CALLER.out.versions)

    /*
     * Taxonomic classification of the contigs with CATPACK
     * This outside of the taxonomical_annotation suboworkflow because it has a dependency with the
     * CGC predicted proteins
    */
    ASSEMBLY_QC.out.assembly_qc_pass
        .join( COMBINED_GENE_CALLER.out.faa )
        .multiMap { meta, contigs, faa ->
            contigs: [meta, contigs]
            proteins: [meta, faa]
        }.set {
            ch_contigs_taxonomic_classification
        }

    CONTIGS_TAXONOMIC_CLASSIFICATION(
        ch_contigs_taxonomic_classification.contigs,
        ch_contigs_taxonomic_classification.proteins,
        [[id: "cat_diamond_db"], file(params.cat_diamond_database, checkIfExists: true)],
        [[id: "cat_taxonomy_db"], file(params.cat_taxonomy_database, checkIfExists: true)],
        params.cat_database_version,
    )
    ch_versions = ch_versions.mix(CONTIGS_TAXONOMIC_CLASSIFICATION.out.versions)

    /*********************/
    /* PROTEINS CHUNKING */
    /*********************/

    // Chunk the fasta into files with at most params.proteins_chunksize sequences
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
        COMBINED_GENE_CALLER.out.faa.join(COMBINED_GENE_CALLER.out.gff),
        ch_protein_chunks,
    )
    ch_versions = ch_versions.mix(FUNCTIONAL_ANNOTATION.out.versions)

    /*
    * Pathway and systems annotations
    */
    PATHWAYS_AND_SYSTEMS(
        ch_protein_chunks,
        ASSEMBLY_QC.out.assembly_qc_pass.join(
            COMBINED_GENE_CALLER.out.faa
        ).join(
            COMBINED_GENE_CALLER.out.gff
        ).join(
            FUNCTIONAL_ANNOTATION.out.interproscan_tsv
        ),
        // DRAM //
        COMBINED_GENE_CALLER.out.faa,
        FUNCTIONAL_ANNOTATION.out.interproscan_tsv,
        FUNCTIONAL_ANNOTATION.out.kegg_orthologs_per_contig_tsv,
        FUNCTIONAL_ANNOTATION.out.dbcan_overview
    )
    ch_versions = ch_versions.mix(PATHWAYS_AND_SYSTEMS.out.versions)

    // Generate giant GFF summary file //
    GFF_SUMMARY(COMBINED_GENE_CALLER.out.gff.join(
        FUNCTIONAL_ANNOTATION.out.interproscan_tsv
        ).join(
            FUNCTIONAL_ANNOTATION.out.eggnog_annotations
        ).join(
            FUNCTIONAL_ANNOTATION.out.dbcan_overview
        ).join(
            FUNCTIONAL_ANNOTATION.out.dbcan_hmm
        ).join(
            PATHWAYS_AND_SYSTEMS.out.sanntis_gff
        ).join(
            PATHWAYS_AND_SYSTEMS.out.antismash_gff
        )
    )
    ch_versions = ch_versions.mix(GFF_SUMMARY.out.versions)

    // Run the genometools GFF validation
    // If there are any errors with the GFF, these will be added to a log file for further inspection,
    // but an invalid GFF won't make the pipeline fail. Doing so would be a waste of compute and not
    // something we could take care of at this point of the execution.
    // This is based on the assumption that we will test the pipeline with loads of different kinds of assemblies and by then
    // we should have ironed out any issues.
    // The log file will be used in the automation system to flag the jobs for inspection.
    // But, I'm up to discuss this @mberacochea
    // TODO: Maybe add a parameter to make the pipeline fail on invalid GFF?
    GT_GFF3VALIDATOR(
        GFF_SUMMARY.out.gff_summary
    )
    ch_versions = ch_versions.mix(GT_GFF3VALIDATOR.out.versions)

    // Collect decontamination TSV files for MultiQC per assembly
    // Use join with remainder: true to handle optional decontamination channels
    ch_per_assembly_files_compressed = ASSEMBLY_QC.out.quast_report_tsv
        .join(ASSEMBLY_QC.out.human_contaminated_contigs_tsv, remainder: true)
        .join(ASSEMBLY_QC.out.phix_contaminated_contigs_tsv, remainder: true)
        .join(ASSEMBLY_QC.out.host_contaminated_contigs_tsv, remainder: true)
        .map { meta, quast, human, phix, host ->
            def files = [quast]
            if (human) files.add(human)
            if (phix) files.add(phix)
            if (host) files.add(host)
            [meta, files]
        }

    // Un-compress for multiqc //
    FIND_UNPIGZ(ch_per_assembly_files_compressed)

    ch_versions = ch_versions.mix(FIND_UNPIGZ.out.versions)

    //
    // Collate and save software versions (it also includes the database versions)
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
        FIND_UNPIGZ.out.file_out,
        common_files,
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    MULTIQC_PER_SAMPLESHEET(
        FIND_UNPIGZ.out.file_out.map { _meta, files -> files }.collect().map { files -> [[id:"samplesheet"], files] },
        common_files,
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    /***********************/
    /* Per assembly report */
    /***********************/

    // This needs to be extended with assemblies that don't have annotations,
    // specifically those without IPS annotations or BGC. We need to retrieve
    // examples of these assemblies from the current backlog database.

    GT_GFF3VALIDATOR.out.success_log
        .join( GT_GFF3VALIDATOR.out.error_log, remainder: true )
        .map { meta, _gff_success, gff_validation_error -> {
                def message = "success"
                if ( gff_validation_error ) {
                    message = "invalid_summary_gff"
                }
                return "${meta.id},${message}"
            }
        }
        .collectFile(name: "analysed_assemblies.csv", storeDir: params.outdir, newLine: true, cache: false)

    // QC Failed  assemblies //
    ASSEMBLY_QC.out.qc_failed_assemblies
        .filter { _meta, qc_failed_message -> qc_failed_message }
        .map { meta, qc_failed_message -> {
                return "${meta.id},${qc_failed_message}"
            }
        }
        .collectFile(name: "qc_failed_assemblies.csv", storeDir: params.outdir, newLine: true, cache: false)

    /*************************/
    /* Downtream samplesheet */
    /************************/

    // VIRIfy samplesheet //
    ASSEMBLY_QC.out.assembly_qc_pass.join(
        COMBINED_GENE_CALLER.out.gff
    ).map {
        meta, assembly, gff -> {
            // We need to handle relative paths
            def outdir_file = file(params.outdir)
            def output_full_path = "${outdir_file.getParent()}/${outdir_file.getName()}"
            return "${meta.id},,,${output_full_path}/${assembly.name},${output_full_path}/${gff.name}"
        }
    }.collectFile(
        name: "virify_samplesheet.csv",
        storeDir: "${params.outdir}/downstream_samplesheets/",
        newLine: true,
        cache: false,
        seed: "id,assembly,fastq_1,fastq_2,proteins"
    )

    emit:
    multiqc_report = MULTIQC_PER_SAMPLESHEET.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
