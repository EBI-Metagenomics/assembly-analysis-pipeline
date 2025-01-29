
/* NF-CORE */
include { DIAMOND_BLASTP                          } from '../../modules/nf-core/diamond/blastp/main'
include { DIAMOND_BLASTP as DIAMOND_FOR_RHEACHEBI } from '../../modules/nf-core/diamond/blastp/main'

/* EBI-METAGENOMICS */
include { INTERPROSCAN                               } from '../../modules/ebi-metagenomics/interproscan/main'
include { EGGNOGMAPPER                               } from '../../modules/ebi-metagenomics/eggnogmapper/main'
include { GENOMEPROPERTIES                           } from '../../modules/ebi-metagenomics/genomeproperties/main'
include { MGNIFYPIPELINESTOOLKIT_RHEACHEBIANNOTATION } from '../../modules/ebi-metagenomics/mgnifypipelinestoolkit/rheachebiannotation/main'

include { GOSLIM_SWF                      } from '../../subworkflows/ebi-metagenomics/goslim_swf/main'

/* LOCAL */
include { HMMER_HMMSCAN as HMMSCAN_KOFAMS } from '../../modules/local/hmmer/hmmscan/main'

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
        "txt", // blast like txt output
        "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sphylums skingdoms sscinames"
    )

    ch_versions = ch_versions.mix(DIAMOND_BLASTP.out.versions)

    /*
     * Assign Rhea and CHEMI tags to the proteins.
     * This is a 2 step process
     * 1 - run diamond against a post-processed UniProt90 + Rhea DB
     * 2 - extract the hits using the mgnif toolkit add rhea annotations script
    */
    DIAMOND_FOR_RHEACHEBI(
        ch_predicted_proteins,
        [["id": "uniref90rheachebi"] , file(params.uniref90rheachebi_diamond_database, checkIfExists: true)],
        "txt", // blast like txt output
        "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    )

    MGNIFYPIPELINESTOOLKIT_RHEACHEBIANNOTATION(
        ch_predicted_proteins.join( DIAMOND_FOR_RHEACHEBI.out.tsv ),
        file(params.uniref90_diamond_database, checkIfExists: true)
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
