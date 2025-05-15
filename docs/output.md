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

The `taxonomy` directory contains output files from taxonomic assignment tools such as `CAT_pack` and `cmsearch`, as well as visualizations generated with `Krona`.

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
- **ERZ12345_SSU.fasta.gz**: This `fasta` file contains all of the matching sequences to a particular rRNA marker gene type from running `cmsearch` on the contigs, being in this case the SSU. As described previously, this file name would be different if a different marker gene was matched, e.g. the suffix could be `_LSU` if it were matching to the LSU model, and you could have multiple files depending on how many marker genes were detected.

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
│   │   ├── ERZ12345_emapper_seed_orthologs.tsv.gz
│   │   └── ERZ12345_emapper_annotations.tsv.gz
│   ├── kegg
│   │   ├── ERZ12345_ko_summary.tsv.gz
│   │   └── ERZ12345_ko_summary.tsv.gz.gzi
│   ├── rhea-reactions
│   │   ├── ERZ12345_proteins2rhea.tsv.gz
│   │   └── ERZ12345_proteins2rhea.tsv.gz.gzi
│   └── dbcan
│       ├── ERZ12345_dbcan_cgc.gff.gz
│       ├── ERZ12345_dbcan_standard_out.tsv.gz
│       ├── ERZ12345_dbcan_overview.tsv.gz
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

- **ERZ12345_emapper_seed_orthologs.tsv.gz**: This `tsv` file contains the seed ortholog matches from running `EggNOG-mapper`.
- **ERZ12345_emapper_annotations.tsv.gz**: This `tsv` file contains the the different annotations from running `EggNOG-mapper`.

#### Output files - kegg

This subdirectory contains the output of running `hmmsearch` on KEGG orthologs, and then summarising their counts.

- **ERZ12345_ko_summary.tsv.gz**: This `tsv` file contains the summary counts of the different KOs that appear in the assembly.
- **ERZ12345_ko_summary.tsv.gz.gzi**: This file is an index for the KO count summary file.

#### Output files - rhea-reactions

This subdirectory contains the output of running `DIAMOND` on a database containing `Rhea` reactions, which are then also linked to their different `ChEBI` accessions.

- **ERZ12345_proteins2rhea.tsv.gz**: This `tsv` file contains the different Rhea IDs and ChEBI reactions in the proteins of the assembly.
- **ERZ12345_proteins2rhea.tsv.gz.gzi**: This file is an index for the Proteins2Rhea output.

#### Output files - dbcan

This subdirectory contains the output of running `run_dbCAN` on the assembly proteins, generating annotations for CAZymes and CAZyme Gene Clusters (CGC).

- **ERZ12345_dbcan_cgc.gff.gz**: This `gff` file is annotated with functional genes for CGCFinder and visualization.
- **ERZ12345_dbcan_standard_out.tsv.gz**: This `tsv` file lists all identified CGCs and their components.
- **ERZ12345_dbcan_overview.tsv.gz**: This `tsv` file contains a summary of identified CAZymes.
- **ERZ12345_dbcan_sub_hmm.tsv.gz**: This file `tsv` file contains the detailed HMMER results.
- **ERZ12345_dbcan_substrates.tsv.gz**: This `tsv` file contains the dbCAN sub-HMM results, including substrate specificity.

### pathways-and-systems

The `pathways-and-systems` directory contains five subdirectories of results, one for each pathway/system annotation tool used by the pipeline. The results range from KEGG pathways, to Biosynthetic Gene Clusters (BGCs) by `antiSMASH` and `SanntiS`. Just like the functional annotations, most of these annotations are on a per-protein basis.

```bash
├── qc
├── cds
├── taxonomy
├── functional-annotation
├── pathways-and-systems
│   ├── antismash
│   │   ├── ERZ12345_antismash.gbk.gz
│   │   ├── ERZ12345_antismash.gff.gz
│   │   ├── ERZ12345_merged.json
│   │   └── ERZ12345_antismash_summary.tsv.gz
│   ├── sanntis
│   │   └── ERZ12345_sanntis_concatenated.gff.gz
│   ├── genome-properties
│   │   ├── ERZ12345_gp.json.gz
│   │   ├── ERZ12345_gp.tsv.gz
│   │   └── ERZ12345_gp.txt.gz
│   ├── kegg-modules
│   │   ├── ERZ12345_kegg_modules_per_contigs.tsv.gz
│   │   ├── ERZ12345_kegg_modules_per_contigs.tsv.gz.gzi
│   │   ├── ERZ12345_kegg_modules_summary.tsv.gz
│   │   └── ERZ12345_kegg_modules_summary.tsv.gz.gzi
│   └── dram-distill
│       ├── ERZ12345_dram.tsv.gz
│       ├── ERZ12345_dram.html.gz
│       ├── ERZ12345_genome_stats.tsv.gz
│       └── ERZ12345_metabolism_summary.xlsx.gz
└── gff
```

#### Output files - antismash

This subdirectory contains the outputs of running `AntiSMASH` on the proteins of the assembly, describing the detected BGCs and outputting them in multiple different formats.

- **ERZ12345_antismash.gbk.gz**: This `gbk` file contains the different AntiSMASH annotations in the GenBank format.
- **ERZ12345_antismash.gff.gz**: This `gff` file contains the different AntiSMASH annotations in the GFF3 format.
- **ERZ12345_merged.json**: This `json` file contains the different AntiSMASH annotations in a JSON object file.
- **ERZ12345_antismash_summary.tsv.gz**: This `tsv` file contains the different AntiSMASH labels in the proteins of the assembly.

#### Output files - sanntis

This subdirectory contains the outputs of running `SanntiS` on the proteins of the assembly, also describing the detected BGCs using this new machine learning-based tool.

- **ERZ12345_sanntis.gff.gz**: This `gff` file contains the different SanntiS annotations in the GFF3 format.

#### Output files - genome-properties

This subdirectory contains the outputs of running `Genome Properties` on the proteins of the assembly, describing different sets of protein signatures that exist in the set of proteins in the assembly, in three different formats.

- **ERZ12345_gp.json.gz**: This `json` file contains the different Genome Properties annotations in a JSON object file.
- **ERZ12345_gp.tsv.gz**: This `tsv` file contains the different Genome Properties annotations in a tab-separated file.
- **ERZ12345_gp.txt.gz**: This `txt` file contains the different Genome Proeprties annotations in a simple text file.

#### Output files - kegg-modules

This subdirectory contains the outputs of computing the different modules and pathsways that the set of KEGG KOs make up in the assembly. Completeness measures for different modules are computed based on this presence of KOs.

- **ERZ12345_kegg_modules_per_contigs.tsv.gz**: This `tsv` file contains the different KEGG modules and their completeness, including the contig they originate from.
- **ERZ12345_kegg_modules_per_contigs.tsv.gz.gzi**: This file is an index for the per-contig KEGG modules file.
- **ERZ12345_kegg_modules_summary.tsv.gz**: This `tsv` file contains a summary of the different KEGG modules and their completeness for the whole assembly, not on a per-contig basis
- **ERZ12345_kegg_modules_summary.tsv.gz.gzi**: This file is an index for the summary KEGG modules file.

#### Output files - dram-distill

This subdirectory contains the outputs of running `DRAM-distill` on some of the outputs of the pipeline, generating analysis and summary reports and visualisations for the annotations on a per-assembly basis, into four different files.

- **ERZ12345_dram.tsv.gz**: This `tsv` file contains the product of `DRAM-distill` for the assembly, including functions that were detected.
- **ERZ12345_dram.html.gz**: This `html` file contains a heatmap visualisation of the detected functions by `DRAM-distill`.
- **ERZ12345_genome_stats.tsv.gz**: This `tsv` file contains a summary of the assembly quality, including number of scaffolds.
- **ERZ12345_metabolism_summary.xlsx.gz**: This `xlsx` file is an Excel spreadsheet containing information about annotations that represent common metabolisms.

### gff

The `gff` directory contains a single file that summarises most of the functional annotations generated by the pipeline on a per-protein basis into a single file.

```bash
├── qc
├── cds
├── taxonomy
├── functional-annotation
├── pathways-and-systems
└── gff
    └── ERZ12345_annotation_summary.gff.gz
```

#### Output files

- **ERZ12345_annotation_summarygz**: This `gff` file contains an expansive and large summary of most of the functional annotations each protein has into a single GFF3 format file.

## Per-study output files

The pipeline generated four different per-study output files that aggregate and summarise data, from successful assembly analysis runs to study-wide `MultiQC` reports. These are stored at the root of the study analysis directory.

### Successfully analyed assemblies

The IDs of assemblies that are successfully analysed are aggregated into a top-level file (`analysed_assemblies.csv `), which looks like this::

```
ERZ12345,success
ERZ56789,success
```

### MultiQC

Just like the pipeline generates MultiQC reports on a per-assembly basis, it also generates the same kind of report but on a per-study basis, which are located in the `multiqc/` directory.

### dram-distill

Just like the pipeline generates `DRAM-distil` outputs on a per-assembly basis, it also generates the same kind of outputs on a per-study basis, which are located in the `dram-distill/` directory.
