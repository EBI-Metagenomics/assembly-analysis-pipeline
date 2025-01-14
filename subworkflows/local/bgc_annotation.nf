/* NF-CORE */
include { ANTISMASH_ANTISMASHLITE } from '../../modules/nf-core/antismash/antismashlite/main'

/* EBI-METAGENOMICS */
include { SANNTIS                 } from '../../modules/ebi-metagenomics/sanntis/main'


workflow BGC_ANNOTATION {

    take:
    ch_predicted_proteins // tule (meta, faa, gff, ips_tsv)

    main:

    ch_versions = Channel.empty()

    antismash_ch = ch_predicted_proteins.multiMap { meta, faa, gff, _ips_tsv ->
        faa: [meta, faa]
        gff: gff
    }

    ANTISMASH_ANTISMASHLITE(
        antismash_ch.faa,
        params.antismash_database,
        params.antismash_installdir,
        antismash_ch.gff
    )

    SANNTIS(
        ch_predicted_proteins.map { meta, _faa, _gff, ips_tsv -> [meta, ips_tsv]}
    )

    emit:
    versions = ch_versions
}
