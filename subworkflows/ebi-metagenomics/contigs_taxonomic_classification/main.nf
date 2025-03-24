include { DIAMOND_BLASTP                     } from '../../../modules/ebi-metagenomics/diamond/blastp/main'
include { CATPACK_CONTIGS                    } from '../../../modules/ebi-metagenomics/catpack/contigs/main'
/* Extension, this will be merge into nf-modules */
include { CATPACK_ADDNAMES                   } from '../../../modules/local/catpack/addnames/main'
include { KRONA_TXT_FROM_CAT_CLASSIFICATION  } from '../../../modules/local/krona_txt_from_cat_classification'
include { KRONA_KTIMPORTTEXT                 } from '../../../modules/local/krona/ktimporttext/main'

// TODO: this is an extension to generate a bgzip file (used on the website)
include { TABIX_BGZIP as TABIX_BGZIP_KRONATXT } from '../../../modules/nf-core/tabix/bgzip/main'


workflow CONTIGS_TAXONOMIC_CLASSIFICATION {
    take:
    contigs     // [ val(meta), path(assembly_fasta) ]
    proteins    // [ val(meta), path(proteins_fasta) ]
    cat_db      // [ val(meta), path(catdb_folder)  ]
    taxonomy_db // [ val(meta), path(cattax_folder) ]

    main:

    ch_versions = Channel.empty()

    /*
    * Note:
    * The CAT tool does not use the provided cat_db, as the alignment is performed by the Diamond step.
    * However, this option is mandatory in the CAT source code.
    */

    DIAMOND_BLASTP(
        proteins,
        [[id: "cat-db"], file("${cat_db[1]}/*.dmnd", checkIfExists: true)],
        6, // blast - txt
        []
    )
    ch_versions = ch_versions.mix(DIAMOND_BLASTP.out.versions.first())

    catpack_input_ch = contigs
        .join( proteins )
        .join( DIAMOND_BLASTP.out.txt )
        .multiMap { meta, contigs_fasta, proteins_faa, diamond_txt ->
            contigs: [meta, contigs_fasta]
            proteins: [meta, proteins_faa]
            diamond_txt: [meta, diamond_txt]
        }

    CATPACK_CONTIGS(
        catpack_input_ch.contigs,
        cat_db,
        taxonomy_db,
        catpack_input_ch.proteins,
        catpack_input_ch.diamond_txt
    )
    ch_versions = ch_versions.mix(CATPACK_CONTIGS.out.versions.first())

    CATPACK_ADDNAMES(
        CATPACK_CONTIGS.out.contig2classification,
        taxonomy_db
    )
    ch_versions = ch_versions.mix(CATPACK_ADDNAMES.out.versions.first())

    KRONA_TXT_FROM_CAT_CLASSIFICATION(
        CATPACK_ADDNAMES.out.txt
    )
    ch_versions = ch_versions.mix(KRONA_TXT_FROM_CAT_CLASSIFICATION.out.versions.first())

    TABIX_BGZIP_KRONATXT(
        KRONA_TXT_FROM_CAT_CLASSIFICATION.out.krona_txt
    )
    ch_versions = ch_versions.mix(TABIX_BGZIP_KRONATXT.out.versions.first())

    KRONA_KTIMPORTTEXT(
        KRONA_TXT_FROM_CAT_CLASSIFICATION.out.krona_txt
    )
    ch_versions = ch_versions.mix(KRONA_KTIMPORTTEXT.out.versions.first())

    emit:
    diamond_blast_tsv         = DIAMOND_BLASTP.out.txt
    contig2classification_tsv = CATPACK_CONTIGS.out.contig2classification
    krona_html                = KRONA_KTIMPORTTEXT.out.html
    versions                  = ch_versions
}
