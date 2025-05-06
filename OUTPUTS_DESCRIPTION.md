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

The `qc` directory contains output files related to the quality control steps of the pipeline, from the contig length filtering by `seqkit`, to the quality assessment done by `QUAST`. The structure of the `qc` directory contains three output files and one output directory:

```bash
├── qc
    ├── ERZ12345_filtered_contigs.fasta.gz
    ├── ERZ12345.tsv
    ├── multiqc_report.html
    └── multiqc_data/
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
    ├── ERZ12345_predicted_orf.ffn.gz
    ├── ERZ12345_predicted_cds.faa.gz
    └── ERZ12345_predicted_cds.gff.gz
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

The `taxonomy` directory contains output files summarising the various taxonomic assignment results using tools like `MAPseq` and `Krona`. The pipeline uses four different reference databases to assign taxonomy, and the outputs in `taxonomy-summary` are split based on the different reference databases:

- SILVA
- PR2
- UNITE
- ITSoneDB

However, there is added complexity in the output structure for two reasons:

- SILVA is actually split into two output directories: SILVA-SSU and SILVA-LSU.
- There are two kinds of taxonomic results for both SILVA-SSU and PR2: one for the hit matches and one for ASVs, with the latter being named as DADA2-SILVA and DADA2-PR2.

#### Output files - hit matches
```bash
├── qc
├── sequence-categorisation
├── amplified-region-inference
├── primer-identification
├── asv
└── taxonomy-summary
    ├── SILVA-SSU
    │   ├── ERR4334351.html
    │   ├── ERR4334351_SILVA-SSU.mseq
    │   ├── ERR4334351_SILVA-SSU.tsv
    │   └── ERR4334351_SILVA-SSU.txt
    ├── PR2
    │   ├── ERR4334351.html
    │   ├── ERR4334351_PR2.mseq
    │   ├── ERR4334351_PR2.tsv
    │   └── ERR4334351_PR2.txt
    ├── UNITE
    │   ├── ERR4334351.html
    │   ├── ERR4334351_UNITE.mseq
    │   ├── ERR4334351_UNITE.tsv
    │   └── ERR4334351_UNITE.txt
    └── ITSoneDB
        ├── ERR4334351.html
        ├── ERR4334351_ITSoneDB.mseq
        ├── ERR4334351_ITSoneDB.tsv
        └── ERR4334351_ITSoneDB.txt
```

All of the different possible subdirectories have the same four files. Taking PR2 as an example:
- **ERR4334351_PR2.mseq**: This `mseq` file contains the raw MAPseq output for every `infernal/cmsearch` match, i.e. each match's taxonomic assignment.
- **ERR4334351_PR2.txt**: This `txt` file contains the Krona text input that is used to generate the Krona HTML file. It contains the distribution of the different taxonomic assignments.
- **ERR4334351.html**: This `html` file contains the Krona HTML file that interactively displays the distribution of the different taxonomic assignments.
- **ERR4334351_UNITE.tsv**: This `tsv` file contains the read count of every taxonomic assignment similar to the Krona txt file, but in a different easier-to-parse format.

#### Output files - ASVs
```bash
├── qc
├── sequence-categorisation
├── amplified-region-inference
├── primer-identification
├── asv
└── taxonomy-summary
    ├── DADA2-SILVA
    │   ├── ERR4334351_16S-V3-V4_DADA2-SILVA_asv_krona_counts.txt
    │   ├── ERR4334351_16S-V3-V4.html
    │   └── ERR4334351_DADA2-SILVA.mseq
    └── DADA2-PR2
        ├── ERR4334351_16S-V3-V4_DADA2-PR2_asv_krona_counts.txt
        ├── ERR4334351_16S-V3-V4.html
        └── ERR4334351_DADA2-PR2.mseq
```

The two different subdirectories have the same three categories of files, two of which are dynamic in naming. Using DADA2-PR2 as an example:
- **ERR4334351_DADA2-PR2.mseq**: This contains the raw MAPseq output for every ASV, i.e. each ASV's taxonomic assignment. This file is not dynamic.
- **ERR4334351_16S-V3-V4_DADA2-PR2_asv_krona_counts.txt**: This file contains the Krona text input that is used to generate the Krona HTML file for ASV results. This file is dynamic, as it will generate it on a per amplified region-basis. This means that in cases where there are two amplified regions, you will have three of these files - one for each reference database, and one for the concatenation of the two.
- **ERR4334351_16S-V3-V4.html**: This file contains the Krona HTMl file that interactively displays the distribution of the different taxonomic assignments for ASV results. This file is dynamic in the exact same way as the Krona text input file, i.e. based on the inferred amplified region(s).

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

| Exclusion message 	|                                                                                         Description                                                                                        	|
|:-----------------:	|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:	|
| `seqfu_fail`        	| Run had an error after running `seqfu check`. Check the log file in `qc/${id}_seqfu.tsv` for the exact reason                                                                             	|
| `sfxhd_fail`        	| Run had an error related to the suffix of the file `_1/_2` not matching the headers inside the fastq file. Check the log file in `qc/${id}_suffix_header_err.json` for the reads at fault 	|
| `libstrat_fail`     	| Run was predicted to likely not be of AMPLICON sequencing based on base-conservation patterns at the beginning of reads                                                                   	|
| `no_reads`          	| Run had no reads left after running `fastp`                                                                                                                                               	|

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

