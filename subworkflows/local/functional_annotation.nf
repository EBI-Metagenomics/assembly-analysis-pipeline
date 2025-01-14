
/* NF-CORE */
include { INTERPROSCAN            } from '../../modules/nf-core/interproscan/main'
// include { DIAMOND_BLASTP       } from '../../modules/nf-core/diamond/blastp/main'

/* EBI-METAGENOMICS */
include { EGGNOGMAPPER            } from '../../modules/ebi-metagenomics/eggnogmapper/main'
include { GENOMEPROPERTIES        } from '../../modules/ebi-metagenomics/genomeproperties/main'


workflow FUNCTIONAL_ANNOTATION {

    take:
    ch_predicted_proteins // tule (meta, faa)

    main:

    ch_versions = Channel.empty()

    // TODO: add chunking

    INTERPROSCAN(
        ch_predicted_proteins,
        params.interproscan_database
    )

    ch_versions = ch_versions.mix(INTERPROSCAN.out.versions)

    EGGNOGMAPPER(
        ch_predicted_proteins,
        [], // tuple val(meta2), path(annotation_hit_table)
        params.eggnog_data_dir,
        params.eggnog_database,
        params.eggnog_diamond_db,
    )

    ch_versions = ch_versions.mix(EGGNOGMAPPER.out.versions)

    GENOMEPROPERTIES(
        INTERPROSCAN.out.tsv
    )

    ch_versions = ch_versions.mix(GENOMEPROPERTIES.out.versions)

    emit:
    interproscan_tsv  = INTERPROSCAN.out.tsv
    // TODO: review this, he gff3 is used with antiSMASH
    interproscan_gff3 = INTERPROSCAN.out.gff3
    versions          = ch_versions
}
