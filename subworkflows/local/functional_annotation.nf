/* NF-CORE */
include { SEQKIT_SPLIT2                           } from '../../modules/nf-core/seqkit/split2/main'
include { DIAMOND_RHEACHEBI                       } from '../../modules/local/diamond_rheachebi'
include { CAT_CAT as CONCATENATE_INTERPROSCAN_TSV } from '../../modules/nf-core/cat/cat/main'

/* EBI-METAGENOMICS */
include { INTERPROSCAN                            } from '../../modules/ebi-metagenomics/interproscan/main'
include { DIAMOND_BLASTP                          } from '../../modules/ebi-metagenomics/diamond/blastp/main'
include { EGGNOGMAPPER                            } from '../../modules/ebi-metagenomics/eggnogmapper/main'
include { GENOMEPROPERTIES                        } from '../../modules/ebi-metagenomics/genomeproperties/main'
include { DBCAN                                   } from '../../modules/ebi-metagenomics/dbcan/dbcan/main'

include { GOSLIM_SWF                              } from '../../subworkflows/ebi-metagenomics/goslim_swf/main'

/* LOCAL */
include { CONCATENATE_INTERPROSCAN_GFFS           } from '../../modules/local/concatenate_interproscan_gffs'
include { EXTRACT_PFAM_COUNTS                     } from '../../modules/local/extract_pfam_counts'
include { HMMER_HMMSCAN as HMMSCAN_KOFAMS         } from '../../modules/local/hmmer/hmmscan/main'
include { KEGGPATHWAYSCOMPLETENESS                } from '../../modules/ebi-metagenomics/keggpathwayscompleteness/main'


workflow FUNCTIONAL_ANNOTATION {
    take:
    ch_predicted_proteins // tule (meta, faa, gff)

    main:

    ch_versions = Channel.empty()

    ch_proteins_faa = ch_predicted_proteins.map { meta, faa, _gff -> [meta, faa] }
    ch_proteins_gff = ch_predicted_proteins.map { meta, _faa, gff -> [meta, gff] }

    // Chunk the fasta into files with at most >= params - TODO: this needs to be a param
    SEQKIT_SPLIT2(
        ch_proteins_faa,
        params.proteins_chunksize,
    )

    ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions)

    def chunked_proteins = SEQKIT_SPLIT2.out.assembly.transpose()

    INTERPROSCAN(
        ch_proteins_faa,
        [file(params.interproscan_database), params.interproscan_database_version],
    )

    ch_versions = ch_versions.mix(INTERPROSCAN.out.versions)

    EGGNOGMAPPER(
        ch_proteins_faa,
        [[], []],
        params.eggnog_data_dir,
        params.eggnog_database,
        params.eggnog_diamond_database,
    )

    ch_versions = ch_versions.mix(EGGNOGMAPPER.out.versions)

    /*
     * We have to concatenate the IPS tsv and gff files for the downstream tools.
     * For the tsv files we can use `cat`, the tsv doesn't have a header row.
     * The GFF it's a bit tricker, that is why we have a custom script to do so.
    */
    CONCATENATE_INTERPROSCAN_TSV(
        INTERPROSCAN.out.tsv.groupTuple()
    )
    ch_versions = ch_versions.mix(CONCATENATE_INTERPROSCAN_TSV.out.versions)

    /*
     * Process the interproscan TSV and extract the Pfam entries counts
    */
    EXTRACT_PFAM_COUNTS(
        CONCATENATE_INTERPROSCAN_TSV.out.file_out
    )
    ch_versions = ch_versions.mix(EXTRACT_PFAM_COUNTS.out.versions)

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
    // TODO: should we chunk?
    DIAMOND_RHEACHEBI(
        ch_proteins_faa,
        file(params.uniref90rhea_diamond_database, checkIfExists: true),
        file(params.rheachebi_mapping_tsv, checkIfExists: true),
        "qseqid sseqid evalue bitscore stitle",
    )

    /*
     * KEGG Orthologous annotation. This step uses hmmscan to annotation the sequences aginst the kofam HMM models.
     * These HMM models have been extended as described -> TODO: link to the mgnify_pipelines_reference_databases pipeline
    */
    // TODO: review the v5 models - they have a description and a GA value
    HMMSCAN_KOFAMS(
        ch_proteins_faa.map { meta, faa ->
            {
                [meta, file("${params.kofam_hmm_database}/*.hmm*", checkIfExists: true), faa, true, true, true]
            }
        }
    )
    ch_versions = ch_versions.mix(HMMSCAN_KOFAMS.out.versions)

    KEGGPATHWAYSCOMPLETENESS(
        ch_proteins_faa.join(HMMSCAN_KOFAMS.out.domain_summary)
    )
    ch_versions = ch_versions.mix(KEGGPATHWAYSCOMPLETENESS.out.versions)

    // TODO: Do we need to consolidate the GFF before DBCan?
    ch_dbcan = ch_proteins_faa
        .join(ch_proteins_gff)
        .multiMap { meta, faa, gff ->
            faa: [meta, faa]
            gff: [meta, gff]
        }

    // TODO: should we chunk?
    DBCAN(
        ch_dbcan.faa,
        ch_dbcan.gff,
        file(params.dbcan_database, checkIfExists: true),
        "protein" // the DBCAN mode
    )
    ch_versions = ch_versions.mix(DBCAN.out.versions)

    emit:
    interproscan_tsv  = CONCATENATE_INTERPROSCAN_TSV.out.file_out
    interproscan_gff3 = CONCATENATE_INTERPROSCAN_GFFS.out.concatenated_gff
    versions          = ch_versions
}
