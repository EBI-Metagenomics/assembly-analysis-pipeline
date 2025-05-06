# MGnify v6 assembly-analysis-pipeline output documentation

This document describes the different outputs produced by the [MGnify v6 assembly-analysis-pipeline](https://github.com/EBI-Metagenomics/assembly-analysis-pipeline), which generates annotation outputs about taxonomy, function, and pathway information for input assemblies.

The pipeline functions on a per-assembly basis, and the outputs are shaped in a similar way. While most of the outputs are on a per-assembly basis, the pipeline does generate some summaries into different top-level summary files which are aggregated from all of the input assemblies. Assuming all the assemblies are from the same study, these summary files can be considered study-level summaries. In this documentation, the different contents of both per-assembly and per-study outputs will be described.

## Per-run output files

There are six general categories of results, which are separated into six different output directories by the pipeline, and each successful assembly should have all six of these directories:

```bash
├── qc
├── cds
├── taxonomy
├── functional-annotation
├── pathways-and-systems
└── gff
```

As these results are per-assembly, most of the outputs use an assembly ID as a prefix. For the purposes of this documentation, we are using the run ID `ERZ12345`.

### qc

The `qc` directory contains output files related to the quality control steps of the pipeline, from the contig length filtering by `seqkit`, to the quality assessment done by `QUAST`. The structure of the `qc` directory contains three output files and one output directory.

```bash
├── qc
│   ├── ERZ12345_filtered_contigs.fasta.gz
│   ├── ERZ12345.tsv
│   ├── multiqc_report.html
│   └── multiqc_data/
├── cds
├── taxonomy
├── functional-annotation
├── pathways-and-systems
└── gff
```

#### Output files

- **ERZ12345_filtered_contigs.fasta.gz**: This `fasta` file contains the filtered contigs after the removal of those that are shorter than 500 bases, and which have a proportion of ambiguous bases higher than 10%.
- **ERZ12345.tsv**: This `tsv` file contains the QUAST summary output, giving an assessment of the quality of the contigs of this assembly.
- **multiqc_report.html**: This `html` file contains the `MultiQC` report for that assembly. It will combine outputs from `QUAST` and the different software versions used by the pipeline.
- **multiqc_data/**: This `directory` contains the input files used by MultiQC to generate its report.

### cds

The `cds` directory contains output files related to the combined gene caller subworkflow, which calls genes using both `Pyrodigal` and `FragGeneScanRs`, which are Python and Rust interfaces to `Prodigal` and `FragGeneScan` respectively. These called genes are then merged, the outputs of which are stored in this directory in three different formats.

```bash
├── qc
├── cds
│   ├── ERZ12345_predicted_orf.ffn.gz
│   ├── ERZ12345_predicted_cds.faa.gz
│   └── ERZ12345_predicted_cds.gff.gz
├── taxonomy
├── functional-annotation
├── pathways-and-systems
└── gff
```

#### Output files

- **ERZ12345_predicted_orf.ffn.gz**: This `fasta` file contains the nucleotide sequences of each predicted Open Reading Frame (ORF) by the combined gene caller, making up the collection of genes detected in the assembly.
- **ERZ12345_predicted_cds.faa.gz**: This `fasta` file contains the amino acid sequences of each predicted ORF, making up the collection of proteins detected in the assembly.
- **ERZ12345_predicted_cds.gff.gz**: This `gff` file contains the different Coding DNA Sequences (CDS) regions in the GFF3 format.

### taxonomy

The `taxonomy` directory contains output files summarising the various taxonomic assignment results using tools like `CAT_pack`, `cmsearch`, and `Krona`. There are four output files in this directory.

```bash
├── qc
├── cds
├── taxonomy
│   ├── ERZ12345_contigs_taxonomy.tsv.gz
│   ├── ERZ12345.krona.txt
│   ├── ERZ12345.html
│   └── ERZ12345_SSU.fasta.gz
├── functional-annotation
├── pathways-and-systems
└── gff
```

#### Output files

- **ERZ12345_contigs_taxonomy.tsv.gz**: This `tsv` file contains the output from `CAT_pack` that describes taxonomy assignments to contigs in the assembly.
- **ERZ12345.krona.txt**: This `txt` file contains the Krona text input that is used to generate the Krona HTML file. It contains the distribution of the different taxonomic assignments from the `CAT_pack` output.
- **ERZ12345.html**: This `html` file contains the Krona HTML file that interactively displays the distribution of the different taxonomic assignments from `CAT_pack`.
- **ERZ12345_SSU.fasta.gz**: This `fasta` file contains all of the matching sequences to a particular rRNA amplicon type from running `cmsearch` on the contigs, being in this case the SSU. As described previously, this file name would be different if a different amplicon was matched, e.g. the suffix could be `_LSU` if it were matching to the LSU model, and you could have multiple files depending on how many marker genes were detected.

### functional-annotation

The `functional-annotation` directory contains seven subdirectories of results, one for each functional annotation tool used in the pipeline. Specifically, the functional annotations in this directory are assigned on a per-protein basis, and range from tools like `InterProScan` to `dbCAN`.

```bash
├── qc
├── cds
├── taxonomy
├── functional-annotation
│   ├── interpro
│   │   ├── ERZ12345_interproscan.tsv.gz
│   │   ├── ERZ12345_interpro_summary.tsv.gz
│   │   └── ERZ12345_interpro_summary.tsv.gz.gzi
│   ├── pfam
│   │   ├── ERZ12345_pfam_summary.tsv.gz
│   │   └── ERZ12345_pfam_summary.tsv.gz.gzi
│   ├── go
│   │   ├── ERZ12345_go_summary.tsv.gz
│   │   ├── ERZ12345_go_summary.tsv.gz.gzi
│   │   ├── ERZ12345_goslim_summary.tsv.gz
│   │   └── ERZ12345_goslim_summary.tsv.gz.gzi
│   ├── eggnog
│   │   ├── ERZ12345_emapper_annotations.tsv.gz
│   │   └── ERZ12345_emapper_seed_orthologs.tsv.gz
│   ├── kegg
│   │   ├── ERZ12345_ko_summary.tsv.gz
│   │   └── ERZ12345_ko_summary.tsv.gz.gzi
│   ├── rhea-reactions
│   │   ├── ERZ12345_proteins2rhea.tsv.gz
│   │   └── ERZ12345_proteins2rhea.tsv.gz.gzi
│   └── dbcan
│       ├── ERZ12345_dbcan_cgc.gff.gz
│       ├── ERZ12345_dbcan_overview.tsv.gz
│       ├── ERZ12345_dbcan_standard_out.tsv.gz
│       ├── ERZ12345_dbcan_sub_hmm.tsv.gz
│       └── ERZ12345_dbcan_substrates.tsv.gz
├── pathways-and-systems
└── gff
```

#### Output files - interpro

This subdirectory contains the output of running InterProScan on the proteins of the assembly.

- **ERZ12345_interproscan.tsv.gz**: This `tsv` file contains the output of InterProScan, containing the different InterPro annotations assigned to the set of proteins in the assembly.
- **ERZ12345_interpro_summary.tsv.gz**: This `tsv` file contains the summary counts of the different InterPro signatures that appear in the assembly.
- **ERZ12345_interpro_summary.tsv.gz.gzi**: This file is an index for the InterPro count summary file.

#### Output files - pfam

This subdirectory contains the output of extracting the Pfam signatures from the InterProScan results.

- **ERZ12345_pfam_summary.tsv.gz**: This `tsv` file contains the summary counts of the different Pfam signatures that appear in the assembly.
- **ERZ12345_pfam_summary.tsv.gz.gzi**: This file is an index for the Pfam count summary file.

#### Output files - go

This subdirectory contains the output of extracting the Gene Ontology (GO) signatures from the InterProScan results. Also, the existing GO terms are summarised using the metagenomics GO-slim ([see here](https://geneontology.org/docs/go-subset-guide/)) for a more domain-specific and concise set of GO terms.

- **ERZ12345_go_summary.tsv.gz**: This `tsv` file contains the summary counts of the different GO signatures that appear in the assembly.
- **ERZ12345_go_summary.tsv.gz.gzi**: This file is an index for the GO count summary file.
- **ERZ12345_goslim_summary.tsv.gz**: This `tsv` file contains the summary counts of the different GO-slim signatures that appear in the assembly.
- **ERZ12345_goslim_summary.tsv.gz.gzi**: This file is an index for the GO-slim count summary file.

#### Output files - eggnog

This subdirectory contains the output of running `EggNOG-mapper` on the proteins of the assembly.

- **ERZ12345_emapper_annotations.tsv.gz**: This `tsv` file contains the summary counts of the different Pfam signatures that appear in the assembly.
- **ERZ12345_emapper_seed_orthologs.tsv.gz**: This file is an index for the Pfam count summary file.

## Per-study output files

The pipeline generated four different per-study output files that aggregate and summarise data, from failed runs to primer validation metadata.

### MultiQC

The pipeline generates two [MultiQC](https://seqera.io/multiqc/) reports: one per-study (`study_multiqc_report.html`), and one per-run (`qc/${id}_multiqc_report.html`). These reports aggregate a few QC statistics from some of the tools run by the pipeline, including:

- fastp
- cutadapt
- DADA2 (as a custom report)

### QC failed runs

The pipeline runs a couple of sanity and QC checks on every input run. In the case where a run fails, it will be added to a top-level report (`qc_failed_runs.csv`) that aggregates the IDs of any other run that failed, along with the particular reason it failed. For example:

```
ERR6093685,no_reads
ERRSFXHDFAIL,sfxhd_fail
ERRSEQFUFAIL,seqfu_fail
SRRLIBSTRATFAIL,libstrat_fail
```

The different exclusion messages are:

| Exclusion message |                                                                                        Description                                                                                        |
| :---------------: | :---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: |
|   `seqfu_fail`    |                                       Run had an error after running `seqfu check`. Check the log file in `qc/${id}_seqfu.tsv` for the exact reason                                       |
|   `sfxhd_fail`    | Run had an error related to the suffix of the file `_1/_2` not matching the headers inside the fastq file. Check the log file in `qc/${id}_suffix_header_err.json` for the reads at fault |
|  `libstrat_fail`  |                                  Run was predicted to likely not be of AMPLICON sequencing based on base-conservation patterns at the beginning of reads                                  |
|    `no_reads`     |                                                                        Run had no reads left after running `fastp`                                                                        |

### QC passed runs

Similarly to runs that fail QC, runs that pass QC are guaranteed to generate results. The IDs of such runs is aggregated into a top-level file (`qc_passed_runs.csv`). For example:

```
SRR17062740,all_results
ERR4334351,all_results
ERRNOASVS,no_asvs
```

An important thing to note is that while a run might succeed at generating results for the closed-reference based method, it might fail at some extra QC checks required for generating results using the ASV method. For this reason, there are two statuses a passed run can have:

- `all_results` - if results for both methods could be generated
- `no_asvs` - if ASV results could not be generated

### Primer validation summary

The pipeline performs inferrence of primer presence and sequence using [PIMENTO](https://github.com/EBI-Metagenomics/PIMENTO/tree/dev). For any runs where a primer was detected, metadata about it will be aggregated into a top-level primer validation summary file (`primer_validation_summary.json`), including its sequence, region, and identification strategy. For example:

```json
[
  {
    "id": "SRR17062740",
    "primers": [
      {
        "name": "F_auto",
        "region": "V4",
        "strand": "fwd",
        "sequence": "ATTCCAGCTCCAATAG",
        "identification_strategy": "auto"
      },
      {
        "name": "R_auto",
        "region": "V4",
        "strand": "rev",
        "sequence": "GACTACGATGGTATNTAATC",
        "identification_strategy": "auto"
      }
    ]
  },
  {
    "id": "ERR4334351",
    "primers": [
      {
        "name": "341F",
        "region": "V3",
        "strand": "fwd",
        "sequence": "CCTACGGGNGGCWGCAG",
        "identification_strategy": "std"
      },
      {
        "name": "805R",
        "region": "V4",
        "strand": "rev",
        "sequence": "GACTACHVGGGTATCTAATCC",
        "identification_strategy": "std"
      }
    ]
  }
]
```

The value of the `identification_strategy` key can either be:

- `std` - Meaning the primer was matched to one of the standard library primers (more reliable)
- `auto` - Meaning the primer was automatically predicted (less reliable)
