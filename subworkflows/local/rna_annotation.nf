/* NF-CORE */
include { SEQKIT_SPLIT2                             } from '../../modules/nf-core/seqkit/split2/main'
include { CAT_CAT as CONCATENATE_CMSEARCH_DEOVERLAP } from '../../modules/nf-core/cat/cat/main'
include { CAT_CAT as CONCATENATE_EASEL_FASTA        } from '../../modules/nf-core/cat/cat/main'

/* EBI-METAGENOMICS */
include { DETECT_RNA      } from '../../subworkflows/ebi-metagenomics/detect_rna/main'
include { EASEL_ESLSFETCH } from '../../modules/ebi-metagenomics/easel/eslsfetch/main'
include { EXTRACTCOORDS   } from '../../modules/ebi-metagenomics/extractcoords/main'

workflow RNA_ANNOTATION {

    take:
    ch_contigs // tuple (meta, fasta)

    main:

    ch_versions = Channel.empty()

    // Chunk the fasta into files with at most params.proteins_chunksize sequences
    SEQKIT_SPLIT2(
        ch_contigs,
        params.bgc_contigs_chunksize // Define a chunk size for this one
    )
    ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions)

    DETECT_RNA(
        SEQKIT_SPLIT2.out.assembly.transpose(),
        file(params.rfam_cm, checkIfExists: true),
        file(params.rfam_claninfo, checkIfExists: true),
        "cmsearch"
    )
    ch_versions = ch_versions.mix(DETECT_RNA.out.versions)

    CONCATENATE_CMSEARCH_DEOVERLAP(
        DETECT_RNA.out.cmsearch_deoverlap_out.groupTuple()
    )
    ch_versions = ch_versions.mix(CONCATENATE_CMSEARCH_DEOVERLAP.out.versions.first())

    EASEL_ESLSFETCH(
        ch_contigs.join( CONCATENATE_CMSEARCH_DEOVERLAP.out.file_out )
    )
    ch_versions = ch_versions.mix(EASEL_ESLSFETCH.out.versions.first())

    CONCATENATE_CMSEARCH_DEOVERLAP.out.file_out.join( EASEL_ESLSFETCH.out.easel_coords ).multiMap { meta, cmsearch_deoverlap_concatenated_out, easel_matches_fasta ->
        cmsearch_deoverlap: [meta, cmsearch_deoverlap_concatenated_out]
        easel_fasta: [meta, easel_matches_fasta]
    }.set {
        ch_extract_coordinates
    }

    EXTRACTCOORDS(
        ch_extract_coordinates.easel_fasta,
        ch_extract_coordinates.cmsearch_deoverlap
    )
    ch_versions = ch_versions.mix(EXTRACTCOORDS.out.versions.first())

    emit:
    ssu_lsu_coords     = EXTRACTCOORDS.out.concat_ssu_lsu_coords
    ssu_fasta          = EXTRACTCOORDS.out.ssu_fasta
    lsu_fasta          = EXTRACTCOORDS.out.lsu_fasta
    versions           = ch_versions
}
