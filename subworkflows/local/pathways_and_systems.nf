/* NF-CORE */
include { SEQKIT_SEQ as SEQKIT_SEQ_BGC                        } from '../../modules/nf-core/seqkit/seq/main'
include { SEQKIT_SPLIT2                                       } from '../../modules/nf-core/seqkit/split2/main'
include { ANTISMASH_ANTISMASHLITE                             } from '../../modules/nf-core/antismash/antismashlite/main'
include { TABIX_BGZIP as TABIX_BGZIP_KEGGPATHWAYSCOMPLETENESS } from '../../modules/nf-core/tabix/bgzip/main'
/* EBI-METAGENOMICS */
include { SANNTIS                      } from '../../modules/ebi-metagenomics/sanntis/main'
include { GENOMEPROPERTIES             } from '../../modules/ebi-metagenomics/genomeproperties/main'

/* LOCAL */
include { ANTISMASH_JSON_TO_GFF        } from '../../modules/local/antismash_json_to_gff'
include { CONCATENATE_GFFS             } from '../../modules/local/concatenate_gffs'
include { KEGGPATHWAYSCOMPLETENESS     } from '../../modules/ebi-metagenomics/keggpathwayscompleteness/main'


workflow PATHWAYS_AND_SYSTEMS {

    take:
    ch_contigs_and_predicted_proteins // tuple (meta, fasta, faa, gff, ips_tsv)
    ch_kegg_orthologs_summary_tsv     // tuple (meta, kos_summary_tsv)

    main:

    ch_versions = Channel.empty()

    KEGGPATHWAYSCOMPLETENESS(
        ch_kegg_orthologs_summary_tsv
    )
    ch_versions = ch_versions.mix(KEGGPATHWAYSCOMPLETENESS.out.versions)

    TABIX_BGZIP_KEGGPATHWAYSCOMPLETENESS(
        KEGGPATHWAYSCOMPLETENESS.out.kegg_pathways
    )
    ch_versions = ch_versions.mix(TABIX_BGZIP_KEGGPATHWAYSCOMPLETENESS.out.versions)

    GENOMEPROPERTIES(
        ch_contigs_and_predicted_proteins.map { meta, _fasta, _faa, _gff, interpro_tsv -> [meta, interpro_tsv] }
    )
    ch_versions = ch_versions.mix(GENOMEPROPERTIES.out.versions)

    /*
    * For BGC the pipeline uses a different min contig size
    */
    SEQKIT_SEQ_BGC(
        ch_contigs_and_predicted_proteins.map {  meta, fasta, _faa, _gff, _ips_tsv -> [ meta, fasta ] }
    )
    ch_versions = ch_versions.mix(SEQKIT_SEQ_BGC.out.versions)

    // Chunk the fasta into files with at most >= params
    SEQKIT_SPLIT2(
        SEQKIT_SEQ_BGC.out.fastx,
        params.bgc_contigs_chunksize
    )
    ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions)

    def ch_chunked_assembly_fasta = SEQKIT_SPLIT2.out.assembly.transpose()

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

    ANTISMASH_JSON_TO_GFF(
        ANTISMASH_ANTISMASHLITE.out.json_results.groupTuple()
    )
    ch_versions = ch_versions.mix(ANTISMASH_JSON_TO_GFF.out.versions)

    CONCATENATE_GFFS(
        ANTISMASH_JSON_TO_GFF.out.antismash_gff.groupTuple()
    )
    ch_versions = ch_versions.mix(CONCATENATE_GFFS.out.versions)

    // TODO: chunk here too, it's taking a long time
    SANNTIS(
        ch_contigs_and_predicted_proteins.map { meta, _fasta, faa, _gff, ips_tsv -> [meta, ips_tsv, [], faa]}
    )
    ch_versions = ch_versions.mix(SANNTIS.out.versions)

    emit:
    versions = ch_versions
}
