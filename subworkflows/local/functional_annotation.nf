/* NF-CORE */
include { SEQKIT_SPLIT2                                   } from '../../modules/nf-core/seqkit/split2/main'
include { CAT_CAT as CONCATENATE_EGGNOGMAPPER_ORTHOLOGS   } from '../../modules/nf-core/cat/cat/main'
include { CAT_CAT as CONCATENATE_EGGNOGMAPPER_ANNOTATIONS } from '../../modules/nf-core/cat/cat/main'
include { CAT_CAT as CONCATENATE_HMMSEARCH_TBLOUT         } from '../../modules/nf-core/cat/cat/main'
include { CAT_CAT as CONCATENATE_INTERPROSCAN_TSV         } from '../../modules/nf-core/cat/cat/main'
include { CSVTK_CONCAT as CONCATENATE_DBCAN_OVERVIEW      } from '../../modules/nf-core/csvtk/concat/main'
include { CSVTK_CONCAT as CONCATENATE_DBCAN_STANDARD_OUT  } from '../../modules/nf-core/csvtk/concat/main'
include { CSVTK_CONCAT as CONCATENATE_DBCAN_SUBSTRATES    } from '../../modules/nf-core/csvtk/concat/main'
include { CSVTK_CONCAT as CONCATENATE_DBCAN_HMMOUT        } from '../../modules/nf-core/csvtk/concat/main'
include { TABIX_BGZIP as TABIX_BGZIP_GO                   } from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_BGZIP as TABIX_BGZIP_GOSLIM               } from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_BGZIP as TABIX_BGZIP_RHEAANDCHEBI         } from '../../modules/nf-core/tabix/bgzip/main'

/* EBI-METAGENOMICS */
include { INTERPROSCAN                             } from '../../modules/ebi-metagenomics/interproscan/main'
include { DBCAN                                    } from '../../modules/ebi-metagenomics/dbcan/dbcan/main'
include { GOSLIM_SWF                               } from '../../subworkflows/ebi-metagenomics/goslim_swf/main'

/* LOCAL */
include { EGGNOGMAPPER_ORTHOLOGS                            } from '../../modules/local/eggnogmapper_orthologs'
include { EGGNOGMAPPER_ANNOTATIONS                          } from '../../modules/local/eggnogmapper_annotations'
include { SEQKIT_FIX                                        } from '../../modules/local/seqkit_fix'
include { CONCATENATE_GFFS as CONCATENATE_INTERPROSCAN_GFFS } from '../../modules/local/concatenate_gffs'
include { CONCATENATE_GFFS as CONCATENATE_DBCAN_GFFS        } from '../../modules/local/concatenate_gffs'
include { PFAM_SUMMARY                                      } from '../../modules/local/pfam_summary'
include { INTERPRO_SUMMARY                                  } from '../../modules/local/interpro_summary'
include { HMMER_HMMSEARCH as HMMSEARCH_KOFAMS               } from '../../modules/nf-core/hmmer/hmmsearch/main'
include { KEGG_ORTHOLOGS_SUMMARY                            } from '../../modules/local/kegg_orthologs_summary'
include { DIAMOND_RHEACHEBI                                 } from '../../modules/local/diamond_rheachebi'

workflow FUNCTIONAL_ANNOTATION {
    take:
    ch_contigs            // tuple (meta, contigs_fasta)
    ch_predicted_proteins // tule (meta, faa, gff)
    ch_protein_chunked    // tuple (meta, faa_chunk)

    main:

    ch_versions = Channel.empty()

    def ch_proteins_faa = ch_predicted_proteins.map { meta, faa, _gff -> [meta, faa] }
    def ch_proteins_gff = ch_predicted_proteins.map { meta, _faa, gff -> [meta, gff] }

    /*
     * InterProScan results are concatenated (TSV and GFF) and passed on to downstream tools:
     * - GO Slim subworkflow
     * - Genomes properties
     * The concatenated TSV is also used to extract the PFAM and InterPro summaries
    */
    INTERPROSCAN(
        ch_protein_chunked,
        [file(params.interproscan_database, checkIfExists: true), params.interproscan_database_version],
    )
    ch_versions = ch_versions.mix(INTERPROSCAN.out.versions)

    // For this one we use cat/cat as the TSV may have spaces that cause csvtk to complain
    // The tsv is generated directly by InterProScan and for the summary we use csvtk fix-quotes
    // So, I think it's safe to use cat/cat but we should revisit at some point
    // TODO: check we should use csvtk fix-quotes on this one
    CONCATENATE_INTERPROSCAN_TSV(
        INTERPROSCAN.out.tsv.groupTuple(),
    )
    ch_versions = ch_versions.mix(CONCATENATE_INTERPROSCAN_TSV.out.versions)

    /*
     * EggNOG mapper results are concatenated (TSV and GFF) and passed on to downstream tools:
     * - GO Slim subworkflow
     * - Genomes properties
     * The concatenated TSV is also used to extract the PFAM and InterPro summaries
    */
    EGGNOGMAPPER_ORTHOLOGS(
        ch_protein_chunked,
        file(params.eggnog_data_dir, checkIfExists: true),
        file(params.eggnog_database, checkIfExists: true),
        file(params.eggnog_diamond_database, checkIfExists: true),
    )
    ch_versions = ch_versions.mix(EGGNOGMAPPER_ORTHOLOGS.out.versions.first())

    EGGNOGMAPPER_ANNOTATIONS(
        EGGNOGMAPPER_ORTHOLOGS.out.orthologs,
        file(params.eggnog_data_dir, checkIfExists: true),
    )
    ch_versions = ch_versions.mix(EGGNOGMAPPER_ANNOTATIONS.out.versions.first())

    CONCATENATE_EGGNOGMAPPER_ORTHOLOGS(
        EGGNOGMAPPER_ORTHOLOGS.out.orthologs.groupTuple()
    )
    ch_versions = ch_versions.mix(CONCATENATE_EGGNOGMAPPER_ORTHOLOGS.out.versions.first())
    EGGNOGMAPPER_ANNOTATIONS.out.annotations.groupTuple().view()
    CONCATENATE_EGGNOGMAPPER_ANNOTATIONS(
        EGGNOGMAPPER_ANNOTATIONS.out.annotations.groupTuple()
    )
    ch_versions = ch_versions.mix(CONCATENATE_EGGNOGMAPPER_ANNOTATIONS.out.versions.first())

    CONCATENATE_INTERPROSCAN_GFFS(
        INTERPROSCAN.out.gff3.groupTuple()
    )
    ch_versions = ch_versions.mix(CONCATENATE_INTERPROSCAN_GFFS.out.versions)

    /*
     * Get GO term and GO-slim term counts out of an InterProScan .tsv output file
     * We ran GOSLIM once per assembly, hence the groupTuple (the IPS results are one per chunk )
    */
    GOSLIM_SWF(
        CONCATENATE_INTERPROSCAN_TSV.out.file_out,
        file(params.go_obo, checkIfExists: true),
        file(params.goslim_ids, checkIfExists: true),
        file(params.go_banding, checkIfExists: true),
    )
    ch_versions = ch_versions.mix(GOSLIM_SWF.out.versions)

    // TODO: remove this once GOSLIM_SWF produces bgzip/index tsvs
    TABIX_BGZIP_GO(
        GOSLIM_SWF.out.go_summary
    )
    ch_versions = ch_versions.mix(TABIX_BGZIP_GO.out.versions)

    TABIX_BGZIP_GOSLIM(
        GOSLIM_SWF.out.goslim_summary
    )
    ch_versions = ch_versions.mix(TABIX_BGZIP_GOSLIM.out.versions)

    /*
     * Assign Rhea and CHEBI tags to the proteins.
     * This is a 2 step process
     * 1 - run diamond against a post-processed UniProt90 + Rhea DB
     * 2 - extract the hits using the mgnif toolkit add rhea annotations script
    */
    DIAMOND_RHEACHEBI(
        ch_proteins_faa,
        file(params.uniref90rhea_diamond_database, checkIfExists: true),
        file(params.rheachebi_mapping_tsv, checkIfExists: true)
    )
    ch_versions = ch_versions.mix(DIAMOND_RHEACHEBI.out.versions)

    TABIX_BGZIP_RHEAANDCHEBI(
        DIAMOND_RHEACHEBI.out.rhea2proteins_tsv
    )
    ch_versions = ch_versions.mix(TABIX_BGZIP_RHEAANDCHEBI.out.versions)

    ch_proteins_gff.combine(ch_protein_chunked, by: 0).multiMap { meta, gff, faa ->
        faa: [meta, faa]
        gff: [meta, gff]
    }.set {
        ch_dbcan
    }

    DBCAN(
        ch_dbcan.faa,
        ch_dbcan.gff,
        [file(params.dbcan_database, checkIfExists: true), params.dbcan_database_version],
        "protein" // mode
    )
    ch_versions = ch_versions.mix(DBCAN.out.versions)

    CONCATENATE_DBCAN_GFFS(
        DBCAN.out.cgc_gff.groupTuple()
    )
    ch_versions = ch_versions.mix(CONCATENATE_DBCAN_GFFS.out.versions)

    CONCATENATE_DBCAN_OVERVIEW(
        DBCAN.out.overview_output.groupTuple(),
        "tsv",
        "tsv",
        true // compress
    )
    ch_versions = ch_versions.mix(CONCATENATE_DBCAN_OVERVIEW.out.versions)

    CONCATENATE_DBCAN_STANDARD_OUT(
        DBCAN.out.cgc_standard_output.groupTuple(),
        "tsv",
        "tsv",
        true // compress
    )
    ch_versions = ch_versions.mix(CONCATENATE_DBCAN_STANDARD_OUT.out.versions)

    CONCATENATE_DBCAN_SUBSTRATES(
        DBCAN.out.substrate_out.groupTuple(),
        "tsv",
        "tsv",
        true // compress
    )
    ch_versions = ch_versions.mix(CONCATENATE_DBCAN_SUBSTRATES.out.versions)

    /*
    * The dbsub_output is a valid TSV, it contains an extra column with no header
    * it looks like this extra column is a duplicated "Coverage". In order to concatenate the
    * tsv with csvtk (which ensures consistency) we run csvtk fix first to adjust the tsvs
    */
    SEQKIT_FIX(DBCAN.out.dbsub_output)

    ch_versions = ch_versions.mix(SEQKIT_FIX.out.versions.first())

    CONCATENATE_DBCAN_HMMOUT(
        SEQKIT_FIX.out.fixed_tsv.groupTuple(),
        "tsv",
        "tsv",
        true // compress
    )
    ch_versions = ch_versions.mix(CONCATENATE_DBCAN_HMMOUT.out.versions)

    /*
    * KEGG orthologs annotation. This step uses hmmsearch to annotation the sequences aginst the kofam HMM models.
    * These HMM models have been extended as described -> TODO: link to the mgnify_pipelines_reference_databases pipeline
    */
    HMMSEARCH_KOFAMS(
        ch_protein_chunked.map { meta, faa ->
            {
                [meta, file(params.kofam_hmm_database, checkIfExists: true), faa, false, true, true] // boolean flags are: write alignment, tblout and domtbl
            }
        }
    )
    ch_versions = ch_versions.mix(HMMSEARCH_KOFAMS.out.versions)

    /*
    * We concatenate the tblout together, this file is not valid as it has the comments too (the ones that start with #)
    * but because the pipeline processes this file subsequently in this pipeline, this doesn't matter
    */
    CONCATENATE_HMMSEARCH_TBLOUT(
        HMMSEARCH_KOFAMS.out.target_summary.groupTuple()
    )
    ch_versions = ch_versions.mix(CONCATENATE_HMMSEARCH_TBLOUT.out.versions)

    //************************************************//
    //                    Summaries                   //
    //***********************************************//

    PFAM_SUMMARY(
        CONCATENATE_INTERPROSCAN_TSV.out.file_out
    )
    ch_versions = ch_versions.mix(PFAM_SUMMARY.out.versions)

    INTERPRO_SUMMARY(
        CONCATENATE_INTERPROSCAN_TSV.out.file_out
    )
    ch_versions = ch_versions.mix(INTERPRO_SUMMARY.out.versions)

    KEGG_ORTHOLOGS_SUMMARY(
        CONCATENATE_HMMSEARCH_TBLOUT.out.file_out,
        file(params.ko_list_txt, checkIfExists: true)
    )
    ch_versions = ch_versions.mix(KEGG_ORTHOLOGS_SUMMARY.out.versions)

    emit:
    dbcan_overview                 = CONCATENATE_DBCAN_OVERVIEW.out.csv
    dbcan_hmm                      = CONCATENATE_DBCAN_HMMOUT.out.csv
    interproscan_tsv               = CONCATENATE_INTERPROSCAN_TSV.out.file_out
    interproscan_gff3              = CONCATENATE_INTERPROSCAN_GFFS.out.concatenated_gff
    eggnog_annotations             = CONCATENATE_EGGNOGMAPPER_ANNOTATIONS.out.file_out
    kegg_orthologs_per_contig_tsv  = KEGG_ORTHOLOGS_SUMMARY.out.ko_per_contig_tsv
    versions                       = ch_versions
}
