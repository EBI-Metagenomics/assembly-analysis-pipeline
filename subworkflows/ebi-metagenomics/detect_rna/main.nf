// General subwf to detect RNA by using provided RFAM models. It uses cmscan or cmsearch.
// Output: deoverlapped table and chosen fasta file with RNA sequences.

// Use cmscan mode if input fasta file is small and models file is quite big (usecase: mags-catalogues-pipeline)
// Important note: .cm file should be cmpress-ed before execution
// Use cmsearch mode if input fasta is massive and models file contains chosen set of models (usecase: ASA)


include { INFERNAL_CMSEARCH           } from '../../../modules/ebi-metagenomics/infernal/cmsearch/main'
include { INFERNAL_CMSCAN             } from '../../../modules/ebi-metagenomics/infernal/cmscan/main'
include { CONVERTCMSCANTOCMSEARCH     } from '../../../modules/ebi-metagenomics/convertcmscantocmsearch/main'
include { CMSEARCHTBLOUTDEOVERLAP     } from '../../../modules/ebi-metagenomics/cmsearchtbloutdeoverlap/main'
// Because we need this SWF to run on chunked contigs, we tool easel to the calling workflow instead
include { EASEL_ESLSFETCH             } from '../../../modules/ebi-metagenomics/easel/eslsfetch/main'

workflow DETECT_RNA {

    take:
    ch_fasta     // channel: [ val(meta), [ fasta ] ]
    rfam         // folder: rfam for cmsearch/cmscan
    claninfo     // file: claninfo for cmsearchtbloutdeoverlap
    mode         // cmsearch/cmscan

    main:

    ch_versions = Channel.empty()
    cmsearch_ch = Channel.empty()

    if ( mode == 'cmsearch' ) {
        INFERNAL_CMSEARCH(
            ch_fasta,
            rfam
        )
        ch_versions = ch_versions.mix(INFERNAL_CMSEARCH.out.versions.first())
        cmsearch_ch = INFERNAL_CMSEARCH.out.cmsearch_tbl
    }
    else if (mode == 'cmscan') {
       INFERNAL_CMSCAN(
            ch_fasta,
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

    emit:
    cmsearch_deoverlap_out = CMSEARCHTBLOUTDEOVERLAP.out.cmsearch_tblout_deoverlapped   // channel: [ val(meta), [ deoverlapped ] ]
    versions               = ch_versions                                                // channel: [ versions.yml ]
}

