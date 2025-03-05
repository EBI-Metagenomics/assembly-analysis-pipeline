/* NF-CORE */
include { ANTISMASH_ANTISMASHLITE } from '../../modules/nf-core/antismash/antismashlite/main'

/* EBI-METAGENOMICS */
include { SANNTIS                 } from '../../modules/ebi-metagenomics/sanntis/main'


workflow BGC_ANNOTATION {

    take:
    ch_contigs_and_predicted_proteins // tule (meta, fasta, faa, gff, ips_tsv)

    main:

    ch_versions = Channel.empty()

    // TODO: is it safe to chunk for antiSMASH
    // TODO: filter those contigs shorter than 1000pb and chunk

    antismash_ch = ch_contigs_and_predicted_proteins.multiMap { meta, fasta, _faa, gff, _ips_tsv ->
        fasta: [meta, fasta]
        gff: gff
    }

    ANTISMASH_ANTISMASHLITE(
        antismash_ch.fasta,
        file(params.antismash_database, checkIfExists: true),
        antismash_ch.gff
    )

    // SANNTIS(
    //     ch_predicted_proteins.map { meta, _faa, _gff, ips_tsv -> [meta, ips_tsv]}
    // )

    emit:
    versions = ch_versions
}
