/* NF-CORE */
include { SEQKIT_SEQ as SEQKIT_SEQ_BGC } from '../../modules/nf-core/seqkit/seq/main'
include { ANTISMASH_ANTISMASHLITE      } from '../../modules/nf-core/antismash/antismashlite/main'

/* EBI-METAGENOMICS */
include { SANNTIS                      } from '../../modules/ebi-metagenomics/sanntis/main'


workflow BGC_ANNOTATION {

    take:
    ch_contigs_and_predicted_proteins // tule (meta, fasta, faa, gff, ips_tsv)

    main:

    ch_versions = Channel.empty()

    /*
    * For BGC the pipeline uses a different min contig size
    */
    SEQKIT_SEQ_BGC(
        ch_contigs_and_predicted_proteins.map {  meta, fasta, _faa, _gff, _ips_tsv -> [ meta, fasta ] }
    )

    ch_versions = ch_versions.mix(SEQKIT_SEQ_BGC.out.versions)

    def ch_chunked_assembly_fasta = SEQKIT_SEQ_BGC.out.fastx.transpose()

    // TODO: Filter also the GFF
    ch_chunked_assembly_fasta.join(ch_contigs_and_predicted_proteins).multiMap { meta, chunked_fasta, _fasta, _faa, gff, ips_tsv ->
        fasta: [meta, chunked_fasta]
        ips_tsv: [meta, ips_tsv]
        gff: gff
    }.set {
        antismash_ch
    }

    ANTISMASH_ANTISMASHLITE(
        antismash_ch.fasta,
        file(params.antismash_database, checkIfExists: true),
        antismash_ch.gff
    )
    ch_versions = ch_versions.mix(ANTISMASH_ANTISMASHLITE.out.versions)

    SANNTIS(
        ch_contigs_and_predicted_proteins.map { meta, _fasta, faa, _gff, ips_tsv -> [meta, ips_tsv, [], faa]}
    )
    ch_versions = ch_versions.mix(SANNTIS.out.versions)

    emit:
    versions = ch_versions
}
