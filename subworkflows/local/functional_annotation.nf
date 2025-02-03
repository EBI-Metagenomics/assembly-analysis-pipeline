
/* NF-CORE */
include { DIAMOND_BLASTP                  } from '../../modules/nf-core/diamond/blastp/main'
include { DIAMOND_RHEACHEBI               } from '../../modules/local/diamond_rheachebi'

/* EBI-METAGENOMICS */
include { INTERPROSCAN                    } from '../../modules/ebi-metagenomics/interproscan/main'
include { EGGNOGMAPPER                    } from '../../modules/ebi-metagenomics/eggnogmapper/main'
include { GENOMEPROPERTIES                } from '../../modules/ebi-metagenomics/genomeproperties/main'

include { GOSLIM_SWF                      } from '../../subworkflows/ebi-metagenomics/goslim_swf/main'

/* LOCAL */
include { HMMER_HMMSCAN as HMMSCAN_KOFAMS } from '../../modules/local/hmmer/hmmscan/main'
include { KEGGPATHWAYSCOMPLETENESS        } from '../../modules/ebi-metagenomics/keggpathwayscompleteness/main'

workflow FUNCTIONAL_ANNOTATION {

    take:
    ch_predicted_proteins // tule (meta, faa)

    main:

    ch_versions = Channel.empty()

    INTERPROSCAN(
        ch_predicted_proteins,
        [ file(params.interproscan_database), params.interproscan_database_version ]
    )

    ch_versions = ch_versions.mix(INTERPROSCAN.out.versions)

    EGGNOGMAPPER(
        ch_predicted_proteins,
        [[], []], // tuple val(meta2), path(annotation_hit_table)
        params.eggnog_data_dir,
        params.eggnog_database,
        params.eggnog_diamond_database,
    )

    ch_versions = ch_versions.mix(EGGNOGMAPPER.out.versions)

    GENOMEPROPERTIES(
        INTERPROSCAN.out.tsv
    )

    ch_versions = ch_versions.mix(GENOMEPROPERTIES.out.versions)

    /*
    * Perform a DIAMOND search against the UniRef90 database.
    * This step identifies UniRef hits and retrieves the corresponding NCBI taxonomy IDs
    */
    // TODO: this is to be reviewed - https://www.ebi.ac.uk/panda/jira/browse/EMG-6975?src=confmacro
    DIAMOND_BLASTP(
        ch_predicted_proteins,
        [["id": "uniref90"] , file(params.uniref90_diamond_database, checkIfExists: true)],
        6, // blast like txt output
        "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sphylums skingdoms sscinames"
    )

    ch_versions = ch_versions.mix(DIAMOND_BLASTP.out.versions)

    // TODO: collate the taxonomy hits from the diamond output, cols - 'staxids sphylums skingdoms sscinames'

    /*
     * Assign Rhea and CHEMI tags to the proteins.
     * This is a 2 step process
     * 1 - run diamond against a post-processed UniProt90 + Rhea DB
     * 2 - extract the hits using the mgnif toolkit add rhea annotations script
    */
    // TODO: We need a small test db @ochkalova probably has one
    // TODO: This diamond module and the RHEACHEBIANNOTATION will be merged into one
    //       That is because the diamond results are just an intermediary result and to save
    //       storage the rheachebi script reads from the stdin
    DIAMOND_RHEACHEBI(
        ch_predicted_proteins,
        file(params.uniref90rhea_diamond_database, checkIfExists: true),
        file(params.rheachebi_mapping_tsv, checkIfExists: true),
        "qseqid sseqid evalue bitscore stitle"
    )

    /*
     * KEGG Orthologous annotation. This step uses hmmscan to annotation the sequences aginst the kofam HMM models.§
     These HMM models have been extended as described TODO: link to the mgnify_pipelines_reference_databases pipeline
     *
    */
    // TODO: these results needs to be processsed
    HMMSCAN_KOFAMS(
        ch_predicted_proteins.map { meta, proteins -> {
                [meta, file("${params.kofam_hmm_database}/*.hmm*", checkIfExists: true), proteins, true, true, true ]
            }
        }
    )

    ch_versions = ch_versions.mix(HMMSCAN_KOFAMS.out.versions)

    KEGGPATHWAYSCOMPLETENESS(
        ch_predicted_proteins.join( HMMSCAN_KOFAMS.out.domain_summary )
    )

    ch_versions = ch_versions.mix(KEGGPATHWAYSCOMPLETENESS.out.versions)

    /*
     * Get GO term and GO-slim term counts out of an InterProScan .tsv output file
    */
    GOSLIM_SWF(
        INTERPROSCAN.out.tsv,
        file("${params.go_obo}", checkIfExists: true),
        file("${params.goslim_ids}", checkIfExists: true),
        file("${params.go_banding}", checkIfExists: true)
    )

    ch_versions = ch_versions.mix(GOSLIM_SWF.out.versions)

    emit:
    interproscan_tsv  = INTERPROSCAN.out.tsv
    interproscan_gff3 = INTERPROSCAN.out.gff3
    versions          = ch_versions
}
