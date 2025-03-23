// General subwf to detect RNA by using provided RFAM models. It uses cmscan or cmsearch.
// Output: deoverlapped table and chosen fasta file with RNA sequences.

// Use cmscan mode if input fasta file is small and models file is quite big (usecase: mags-catalogues-pipeline)
// Important note: .cm file should be cmpress-ed before execution
// Use cmsearch mode if input fasta is massive and models file contains chosen set of models (usecase: ASA)


include { INFERNAL_CMSEARCH           } from '../../../modules/ebi-metagenomics/infernal/cmsearch/main'
include { INFERNAL_CMSCAN             } from '../../../modules/ebi-metagenomics/infernal/cmscan/main'
include { CONVERTCMSCANTOCMSEARCH     } from '../../../modules/ebi-metagenomics/convertcmscantocmsearch/main'
include { CMSEARCHTBLOUTDEOVERLAP     } from '../../../modules/ebi-metagenomics/cmsearchtbloutdeoverlap/main'
include { EASEL_ESLSFETCH             } from '../../../modules/ebi-metagenomics/easel/eslsfetch/main'

// Extension to join the infernal_cmsearch tables as this pipeline chunks the contigs
include { CAT_CAT                     } from '../../../modules/nf-core/cat/cat/main'

workflow DETECT_RNA {

    take:
    ch_fasta  // channel: [ val(meta), [ fasta_chunks... ], fasta ]
    rfam      // folder: rfam for cmsearch/cmscan
    claninfo  // file: claninfo for cmsearchtbloutdeoverlap
    mode      // cmsearch/cmscan

    main:

    ch_versions = Channel.empty()
    cmsearch_ch = Channel.empty()

    ch_chunked_fasta = ch_fasta.map { meta, _fasta, chunks -> [meta, chunks] }.transpose()
    ch_nonchuncked_fasta = ch_fasta.map { meta, fasta, _chunks -> [meta, fasta] }

    if ( mode == 'cmsearch' ) {
        INFERNAL_CMSEARCH(
            ch_chunked_fasta,
            rfam
        )
        ch_versions = ch_versions.mix(INFERNAL_CMSEARCH.out.versions.first())

        // This is an extension of this SWF for this particular pipeline.
        // For performance and speed purposes, we chunk the FASTA with the contigs.
        // To further parallelize this, it causes issues downstream; different
        // chunks from cmsearch end in easlfecth, causing it to fail. So, we
        // just join them after Infernal to make the code work and to keep it simple.
        CAT_CAT(
            INFERNAL_CMSEARCH.out.cmsearch_tbl.groupTuple()
        )
        cmsearch_ch = CAT_CAT.out.file_out

        ch_versions = ch_versions.mix(CAT_CAT.out.versions.first())
    }
    else if (mode == 'cmscan') {
       INFERNAL_CMSCAN(
            ch_chunked_fasta,
            rfam
       )
       ch_versions = ch_versions.mix(INFERNAL_CMSCAN.out.versions.first())

       CONVERTCMSCANTOCMSEARCH(INFERNAL_CMSCAN.out.cmscan_tbl)
       ch_versions = ch_versions.mix(CONVERTCMSCANTOCMSEARCH.out.versions.first())

       cmsearch_ch = CONVERTCMSCANTOCMSEARCH.out.cmsearch_tblout
    }

    CMSEARCHTBLOUTDEOVERLAP(
        cmsearch_ch,
        claninfo
    )
    ch_versions = ch_versions.mix(CMSEARCHTBLOUTDEOVERLAP.out.versions.first())

    ch_easel = ch_nonchuncked_fasta
                .join(CMSEARCHTBLOUTDEOVERLAP.out.cmsearch_tblout_deoverlapped)

    EASEL_ESLSFETCH(
        ch_easel
    )
    ch_versions = ch_versions.mix(EASEL_ESLSFETCH.out.versions.first())

    emit:
    cmsearch_deoverlap_out = CMSEARCHTBLOUTDEOVERLAP.out.cmsearch_tblout_deoverlapped   // channel: [ val(meta), [ deoverlapped ] ]
    easel_out              = EASEL_ESLSFETCH.out.easel_coords                           // channel: [ val(meta), [ fasta ] ]
    versions               = ch_versions                                                // channel: [ versions.yml ]
}

