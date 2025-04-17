/* NF-CORE */
include { SEQKIT_SEQ as SEQKIT_SEQ_BGC                                   } from '../../modules/nf-core/seqkit/seq/main'
include { SEQKIT_SPLIT2                                                  } from '../../modules/nf-core/seqkit/split2/main'
include { ANTISMASH_ANTISMASHLITE                                        } from '../../modules/nf-core/antismash/antismashlite/main'
include { TABIX_BGZIP as TABIX_BGZIP_KEGGPATHWAYSCOMPLETENESS            } from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_BGZIP as TABIX_BGZIP_KEGGPATHWAYSCOMPLETENESS_PER_CONTIG } from '../../modules/nf-core/tabix/bgzip/main'

/* EBI-METAGENOMICS */
include { SANNTIS                      } from '../../modules/ebi-metagenomics/sanntis/main'
include { GENOMEPROPERTIES             } from '../../modules/ebi-metagenomics/genomeproperties/main'
include { KEGGPATHWAYSCOMPLETENESS     } from '../../modules/ebi-metagenomics/keggpathwayscompleteness/main'

/* LOCAL */
include { ANTISMASH_JSON_TO_GFF                          } from '../../modules/local/antismash_json_to_gff'
include { CONCATENATE_GFFS as CONCATENATE_ANTISMASH_GFFS } from '../../modules/local/concatenate_gffs'
include { CONCATENATE_GFFS as CONCATENATE_SANNTIS_GFFS   } from '../../modules/local/concatenate_gffs'
include { ANTISMASH_SUMMARY                              } from '../../modules/local/antismash_summary'
include { SANNTIS_SUMMARY                                } from '../../modules/local/sanntis_summary'

include { DRAM_DISTILL_SWF                               } from '../../subworkflows/local/dram_distill_swf'

workflow PATHWAYS_AND_SYSTEMS {

    take:
    // Chunked proteins, used in the functional_annotation mostly, we need this for SanntiS
    ch_protein_chunks // tuple (meta, faa_chunk)

    // fasta: contigs
    // faa: CGC predictions faa
    // gff: CGC predictions gff
    // ips_ts: interpsocan concatenated tsv (all the IPS annotations)
    ch_contigs_and_predicted_proteins // tuple (meta, fasta, faa, gff, ips_tsv)


    // Non-chunked proteins
    ch_proteins       // tuple (meta, faa)

    // InterProScan concatenated TSV
    ch_interproscan_tsv

    // KO per contig aggregated - single file per assembly
    ch_kegg_orthologs_per_contig_tsv     // tuple (meta, kos_per_contig_tsv)

    // DBCan concatenated overview
    ch_dbcan_overview                 // tuple (meta, dbcan_overview_tsv)

    main:

    ch_versions = Channel.empty()

    KEGGPATHWAYSCOMPLETENESS(
        ch_kegg_orthologs_per_contig_tsv
    )
    ch_versions = ch_versions.mix(KEGGPATHWAYSCOMPLETENESS.out.versions)

    TABIX_BGZIP_KEGGPATHWAYSCOMPLETENESS(
        KEGGPATHWAYSCOMPLETENESS.out.kegg_pathways
    )
    ch_versions = ch_versions.mix(TABIX_BGZIP_KEGGPATHWAYSCOMPLETENESS.out.versions)

    TABIX_BGZIP_KEGGPATHWAYSCOMPLETENESS_PER_CONTIG(
        KEGGPATHWAYSCOMPLETENESS.out.kegg_pathways_per_contig
    )
    ch_versions = ch_versions.mix(TABIX_BGZIP_KEGGPATHWAYSCOMPLETENESS_PER_CONTIG.out.versions)

    GENOMEPROPERTIES(
        ch_contigs_and_predicted_proteins.map { meta, _fasta, _faa, _gff, interpro_tsv -> [meta, interpro_tsv] }
    )
    ch_versions = ch_versions.mix(GENOMEPROPERTIES.out.versions)

    /*****************************************************************************************
    * For Biosynthetic Gene Clusters (BGC), the pipeline uses a different minimum contig size.
    ******************************************************************************************/
    SEQKIT_SEQ_BGC(
        ch_contigs_and_predicted_proteins.map {  meta, fasta, _faa, _gff, _ips_tsv -> [ meta, fasta ] }
    )
    ch_versions = ch_versions.mix(SEQKIT_SEQ_BGC.out.versions)

    // Chunk the fasta into files with at most params.bgc_contigs_chunksize sequences
    SEQKIT_SPLIT2(
        SEQKIT_SEQ_BGC.out.fastx,
        params.bgc_contigs_chunksize
    )
    ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions)

    /*******************************************************************************************/
    /* Rearrange the channel. We need to create a channel so that                              */
    /* each chunk of the FASTA has the GFF and IPS TSV (these two are for the whole assembly). */
    /*******************************************************************************************/
    def ch_chunked_assembly_fasta = SEQKIT_SPLIT2.out.assembly.transpose()

    // Note: if I do "def bgc_channel = .. combine(...)"" I get this error:
    // def bgc_channel = ch_contigs_and_predicted_proteins.combine(ch_chunked_assembly_fasta, by: 0)
    // weird, that is why there if def ... = channel.empty() and then I assign it
    def antismash_channel = channel.empty()
    antismash_channel = ch_contigs_and_predicted_proteins.combine(ch_chunked_assembly_fasta, by: 0)
        .map { meta, _all_contigs_fasta, _faa, gff, _ips_tsv, contigs_chunk -> [meta, contigs_chunk, gff] }

    ANTISMASH_ANTISMASHLITE(
        antismash_channel,
        file(params.antismash_database, checkIfExists: true)
    )
    ch_versions = ch_versions.mix(ANTISMASH_ANTISMASHLITE.out.versions)

    ANTISMASH_JSON_TO_GFF(
        ANTISMASH_ANTISMASHLITE.out.json_results
    )
    ch_versions = ch_versions.mix(ANTISMASH_JSON_TO_GFF.out.versions.first())

    CONCATENATE_ANTISMASH_GFFS(
        ANTISMASH_JSON_TO_GFF.out.antismash_gff.groupTuple()
    )
    ch_versions = ch_versions.mix(CONCATENATE_ANTISMASH_GFFS.out.versions)

    // TODO: cat gbk and json (?)

    //ANTISMASH_SUMMARY(CONCATENATE_ANTISMASH_GFFS.out.concatenated_gff)
    //ch_versions = ch_versions.mix(ANTISMASH_SUMMARY.out.versions)

    /*************************************************************************************/
    /* Rearrange the channel. We need to create a channel so that                        */
    /* each FAA chunk the IPS TSV (this is the concatenated IPS for the whole assembly). */
    /* We feed in the whole IPS TSV to avoid having to chunk it and match it the protein */
    /* chunks, this wastes some memory as SanntiS loads the whole TSV in memory but this */
    /* should be problematic as the TSV is relativly small (<500MB)                      */
    /*************************************************************************************/
    def sanntis_channel = channel.empty()
    // Note: same weirdness ass antismash_channel
    sanntis_channel = ch_contigs_and_predicted_proteins.combine(ch_protein_chunks, by: 0)
        .map { meta, _all_contigs_fasta, _faa, _gff, ips_tsv, faa_chunk -> [meta, ips_tsv, [], faa_chunk] }

    SANNTIS(
        sanntis_channel
    )
    ch_versions = ch_versions.mix(SANNTIS.out.versions)


    CONCATENATE_SANNTIS_GFFS(
        SANNTIS.out.gff.groupTuple()
    )
    ch_versions = ch_versions.mix(CONCATENATE_SANNTIS_GFFS.out.versions)


    //SANNTIS_SUMMARY(CONCATENATE_SANNTIS_GFFS.out.concatenated_gff)
    //ch_versions = ch_versions.mix(SANNTIS_SUMMARY.out.versions)


    /*
    * DRAM distill - per assembly and for the whole samplesheet
    */
    DRAM_DISTILL_SWF(
        ch_proteins,
        KEGGPATHWAYSCOMPLETENESS.out.kos_aggregated_by_contig, // This is the aggregated ko per contig file
        ch_interproscan_tsv,
        ch_dbcan_overview
    )
    ch_versions = ch_versions.mix(DRAM_DISTILL_SWF.out.versions)

    emit:
    versions = ch_versions
    sanntis_gff = CONCATENATE_SANNTIS_GFFS.out.concatenated_gff
    antismash_gff = CONCATENATE_ANTISMASH_GFFS.out.concatenated_gff
}
