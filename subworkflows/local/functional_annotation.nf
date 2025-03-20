/* NF-CORE */
include { SEQKIT_SPLIT2                                   } from '../../modules/nf-core/seqkit/split2/main'
include { DIAMOND_RHEACHEBI                               } from '../../modules/local/diamond_rheachebi'
include { CAT_CAT as CONCATENATE_EGGNOGMAPPER_ORTHOLOGOUS } from '../../modules/nf-core/cat/cat/main'
include { CAT_CAT as CONCATENATE_EGGNOGMAPPER_ANNOTATIONS } from '../../modules/nf-core/cat/cat/main'
include { CAT_CAT as CONCATENATE_INTERPROSCAN_TSV         } from '../../modules/nf-core/cat/cat/main'
include { CAT_CAT as CONCATENATE_DBCAN_OVERVIEW           } from '../../modules/nf-core/cat/cat/main'
include { CAT_CAT as CONCATENATE_HMMSEARCH_TBLOUT         } from '../../modules/nf-core/cat/cat/main'

/* EBI-METAGENOMICS */
include { INTERPROSCAN                             } from '../../modules/ebi-metagenomics/interproscan/main'
include { EGGNOGMAPPER as EGGNOGMAPPER_ORTHOLOGS   } from '../../modules/ebi-metagenomics/eggnogmapper/main'
include { EGGNOGMAPPER as EGGNOGMAPPER_ANNOTATIONS } from '../../modules/ebi-metagenomics/eggnogmapper/main'
include { GENOMEPROPERTIES                         } from '../../modules/ebi-metagenomics/genomeproperties/main'
include { DBCAN                                    } from '../../modules/ebi-metagenomics/dbcan/dbcan/main'
include { GOSLIM_SWF                               } from '../../subworkflows/ebi-metagenomics/goslim_swf/main'

/* LOCAL */
include { CONCATENATE_GFFS as CONCATENATE_INTERPROSCAN_GFFS } from '../../modules/local/concatenate_gffs'
include { CONCATENATE_GFFS as CONCATENATE_DBCAN_GFFS        } from '../../modules/local/concatenate_gffs'
include { PFAM_SUMMARY                                      } from '../../modules/local/pfam_summary'
include { INTERPRO_SUMMARY                                  } from '../../modules/local/interpro_summary'
include { HMMER_HMMSEARCH as HMMSEARCH_KOFAMS               } from '../../modules/nf-core/hmmer/hmmsearch/main'
include { KEGG_ORTHOLOGS_SUMMARY                            } from '../../modules/local/kegg_orthologs_summary'
include { KEGGPATHWAYSCOMPLETENESS                          } from '../../modules/ebi-metagenomics/keggpathwayscompleteness/main'


workflow FUNCTIONAL_ANNOTATION {
    take:
    ch_predicted_proteins // tule (meta, faa, gff)

    main:

    ch_versions = Channel.empty()

    def ch_proteins_faa = ch_predicted_proteins.map { meta, faa, _gff -> [meta, faa] }
    def ch_proteins_gff = ch_predicted_proteins.map { meta, _faa, gff -> [meta, gff] }

    /*********************/
    /* PROTEINS CHUNKING */
    /*********************/

    // Chunk the fasta into files with at most >= params
    SEQKIT_SPLIT2(
        ch_proteins_faa,
        params.proteins_chunksize,
    )

    ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions)

    def ch_protein_splits = SEQKIT_SPLIT2.out.assembly.transpose()

    /*
     * InterProScan results are concatenated (TSV and GFF) and passed on to downstream tools:
     * - GO Slim subworkflow
     * - Genomes properties
     * The concatenated TSV is also used to extract the PFAM and InterPro summaries
    */
    INTERPROSCAN(
        ch_protein_splits,
        [file(params.interproscan_database, checkIfExists: true), params.interproscan_database_version],
    )
    ch_versions = ch_versions.mix(INTERPROSCAN.out.versions)

    CONCATENATE_INTERPROSCAN_TSV(
        INTERPROSCAN.out.tsv.groupTuple()
    )
    ch_versions = ch_versions.mix(CONCATENATE_INTERPROSCAN_TSV.out.versions)

    /*
     * EggNOG mapper results are concatenated (TSV and GFF) and passed on to downstream tools:
     * - GO Slim subworkflow
     * - Genomes properties
     * The concatenated TSV is also used to extract the PFAM and InterPro summaries
    */
    EGGNOGMAPPER_ORTHOLOGS(
        ch_protein_splits,
        [[], []],
        params.eggnog_data_dir,
        params.eggnog_database,
        params.eggnog_diamond_database,
    )
    ch_versions = ch_versions.mix(EGGNOGMAPPER_ORTHOLOGS.out.versions.first())

    EGGNOGMAPPER_ANNOTATIONS(
        [[], []],
        EGGNOGMAPPER_ORTHOLOGS.out.orthologs,
        params.eggnog_data_dir,
        [],
        []
    )
    ch_versions = ch_versions.mix(EGGNOGMAPPER_ANNOTATIONS.out.versions.first())

    CONCATENATE_EGGNOGMAPPER_ORTHOLOGOUS(
        EGGNOGMAPPER_ORTHOLOGS.out.orthologs.groupTuple()
    )
    ch_versions = ch_versions.mix(CONCATENATE_EGGNOGMAPPER_ORTHOLOGOUS.out.versions.first())


    CONCATENATE_EGGNOGMAPPER_ANNOTATIONS(
        EGGNOGMAPPER_ANNOTATIONS.out.annotations.groupTuple()
    )
    ch_versions = ch_versions.mix(CONCATENATE_EGGNOGMAPPER_ANNOTATIONS.out.versions.first())

    CONCATENATE_INTERPROSCAN_GFFS(
        INTERPROSCAN.out.gff3.groupTuple()
    )
    ch_versions = ch_versions.mix(CONCATENATE_INTERPROSCAN_GFFS.out.versions)

    GENOMEPROPERTIES(
        CONCATENATE_INTERPROSCAN_TSV.out.file_out
    )
    ch_versions = ch_versions.mix(GENOMEPROPERTIES.out.versions)

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

    /*
     * DBCan
    */
    ch_protein_splits.join( ch_proteins_gff ).multiMap { meta, faa, gff ->
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

    // TODO: the header is repeated here
    CONCATENATE_DBCAN_OVERVIEW(
        DBCAN.out.overview_output.groupTuple()
    )
    ch_versions = ch_versions.mix(CONCATENATE_DBCAN_OVERVIEW.out.versions)

    /*
     * KEGG Orthologous annotation. This step uses hmmscan to annotation the sequences aginst the kofam HMM models.
     * These HMM models have been extended as described -> TODO: link to the mgnify_pipelines_reference_databases pipeline
    */
    HMMSEARCH_KOFAMS(
        ch_protein_splits.map { meta, faa ->
            {
                [meta, file(params.kofam_hmm_database, checkIfExists: true), faa, false, true, true] // boolean flags are: write alingment, tblout and domtbl
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
        CONCATENATE_HMMSEARCH_TBLOUT.out.file_out
    )
    ch_versions = ch_versions.mix(KEGG_ORTHOLOGS_SUMMARY.out.versions)

    KEGGPATHWAYSCOMPLETENESS(
        KEGG_ORTHOLOGS_SUMMARY.out.ko_per_contig_tsv
    )
    ch_versions = ch_versions.mix(KEGGPATHWAYSCOMPLETENESS.out.versions)

    emit:
    interproscan_tsv  = CONCATENATE_INTERPROSCAN_TSV.out.file_out
    interproscan_gff3 = CONCATENATE_INTERPROSCAN_GFFS   .out.concatenated_gff
    versions          = ch_versions
}
