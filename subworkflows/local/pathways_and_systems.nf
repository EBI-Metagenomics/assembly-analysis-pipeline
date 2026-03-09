/* NF-CORE */
include { SEQKIT_SEQ as SEQKIT_SEQ_BGC                                   } from '../../modules/nf-core/seqkit/seq/main'
include { SEQKIT_SPLIT2                                                  } from '../../modules/nf-core/seqkit/split2/main'
include { ANTISMASH_ANTISMASH                                            } from '../../modules/nf-core/antismash/antismash/main'
include { TABIX_BGZIP as TABIX_BGZIP_KEGGPATHWAYSCOMPLETENESS            } from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_BGZIP as TABIX_BGZIP_KEGGPATHWAYSCOMPLETENESS_PER_CONTIG } from '../../modules/nf-core/tabix/bgzip/main'
include { CAT_CAT as CONCATENATE_ANTISMASH_GBK                           } from '../../modules/nf-core/cat/cat/main'

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
include { MERGE_ANTISMASH_JSON                           } from '../../modules/local/merge_antismash_json'
include { FILTER_IPS_AND_FAA_BY_CONTIGS                  } from '../../modules/local/filter_ips_and_faa_by_contigs'

include { DRAM_DISTILL_SWF                               } from '../../subworkflows/local/dram_distill_swf'

workflow PATHWAYS_AND_SYSTEMS {

    take:
    // This contains all the proteins (faa and gff) and their IPS annotations
    // fasta: contigs
    // faa: CGC predictions faa
    // gff: CGC predictions gff
    // ips_tsv: interproscan concatenated tsv
    ch_contigs_and_predicted_proteins // tuple (meta, fasta, faa, gff, ips_tsv)

    // KO per contig aggregated - single file per assembly
    ch_kegg_orthologs_per_contig_tsv     // tuple (meta, kos_per_contig_tsv)

    // DBCan concatenated overview
    ch_dbcan_overview                 // tuple (meta, dbcan_overview_tsv)

    main:

    ch_versions = channel.empty()

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

    // Split the contigs FASTA by total nucleotides per chunk
    SEQKIT_SPLIT2(
        SEQKIT_SEQ_BGC.out.fastx,
        params.bgc_analysis_fasta_chunksize,  // length: max number of nucleotides per chunk
        []                                    // size: [disabled] max number of sequences per chunk
    )
    ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions)

    /*******************************************************************************************/
    /* Rearrange the channel. We need to create a channel so that                              */
    /* each chunk of the FASTA has the GFF and IPS TSV (these two are for the whole assembly). */
    /*******************************************************************************************/
    ch_antismash = ch_contigs_and_predicted_proteins
        .combine(SEQKIT_SPLIT2.out.assembly.transpose(), by: 0)
        .map { meta, _all_contigs_fasta, _faa, gff, _ips_tsv, contigs_chunk -> [meta, contigs_chunk, gff] }

    ANTISMASH_ANTISMASH(
        ch_antismash,
        file(params.antismash_database, checkIfExists: true),
        params.antismash_database_version
    )
    ch_versions = ch_versions.mix(ANTISMASH_ANTISMASH.out.versions)

    ANTISMASH_JSON_TO_GFF(
        ANTISMASH_ANTISMASH.out.json_results
    )
    ch_versions = ch_versions.mix(ANTISMASH_JSON_TO_GFF.out.versions.first())

    CONCATENATE_ANTISMASH_GFFS(
        ANTISMASH_JSON_TO_GFF.out.antismash_gff.groupTuple()
    )
    ch_versions = ch_versions.mix(CONCATENATE_ANTISMASH_GFFS.out.versions)

    CONCATENATE_ANTISMASH_GBK(
        ANTISMASH_ANTISMASH.out.gbk_input.groupTuple()
    )
    ch_versions = ch_versions.mix(CONCATENATE_ANTISMASH_GBK.out.versions)

    MERGE_ANTISMASH_JSON(
        ANTISMASH_ANTISMASH.out.json_results.groupTuple()
    )
    ch_versions = ch_versions.mix(MERGE_ANTISMASH_JSON.out.versions)

    // AntiSMASH may not produce any results, the GFF files maybe "empty" .. just the one line (the GFF header)
    // The summary won't be generated in that case
    ANTISMASH_SUMMARY(
        CONCATENATE_ANTISMASH_GFFS.out.concatenated_gff
    )
    ch_versions = ch_versions.mix(ANTISMASH_SUMMARY.out.versions)

    // For each chunk we filter the full IPS TSV and FAA to only the proteins belonging
    // to contigs in that chunk, then run SANNTIS in parallel and merge the results.
    ch_sanntis_chunks = ch_contigs_and_predicted_proteins
        .combine(SEQKIT_SPLIT2.out.assembly.transpose(), by: 0)
        .map { meta, _all_contigs_fasta, faa, _gff, ips_tsv, contigs_chunk ->
            [meta, contigs_chunk, ips_tsv, faa]
        }

    FILTER_IPS_AND_FAA_BY_CONTIGS(ch_sanntis_chunks)
    ch_versions = ch_versions.mix(FILTER_IPS_AND_FAA_BY_CONTIGS.out.versions)

    SANNTIS(
        FILTER_IPS_AND_FAA_BY_CONTIGS.out.filtered
            .map { meta, ips, faa -> [meta, ips, [], faa] }
    )
    ch_versions = ch_versions.mix(SANNTIS.out.versions)

    CONCATENATE_SANNTIS_GFFS(
        SANNTIS.out.gff.groupTuple()
    )
    ch_versions = ch_versions.mix(CONCATENATE_SANNTIS_GFFS.out.versions)

    SANNTIS_SUMMARY(
        CONCATENATE_SANNTIS_GFFS.out.concatenated_gff
    )
    ch_versions = ch_versions.mix(SANNTIS_SUMMARY.out.versions)

    /*
    * DRAM distill - per assembly and for the whole samplesheet
    */
    DRAM_DISTILL_SWF(
        ch_contigs_and_predicted_proteins.map { meta, _fasta, faa, _gff, _ips -> [meta, faa]},
        KEGGPATHWAYSCOMPLETENESS.out.kos_aggregated_by_contig, // This is the aggregated ko per contig file
        ch_contigs_and_predicted_proteins.map { meta, _fasta, _faa, _gff, ips -> [meta, ips]},
        ch_dbcan_overview
    )
    ch_versions = ch_versions.mix(DRAM_DISTILL_SWF.out.versions)

    emit:
    versions = ch_versions
    sanntis_gff = CONCATENATE_SANNTIS_GFFS.out.concatenated_gff
    antismash_gff = CONCATENATE_ANTISMASH_GFFS.out.concatenated_gff
}
