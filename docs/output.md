# MGnify v6 assembly-analysis-pipeline output documentation

This document describes the different outputs produced by the [MGnify v6 assembly-analysis-pipeline](https://github.com/EBI-Metagenomics/assembly-analysis-pipeline), which annotates input assemblies with taxonomy, function, and pathway information.

The pipeline operates on a per-assembly basis, and its outputs generally follow the same structure. While most outputs are generated individually for each assembly, the pipeline also produces several summary files that aggregate data across all input assemblies. If all assemblies belong to the same study, these summary files can be considered study-level outputs. This documentation outlines the contents of both the per-assembly and per-study outputs.

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

As these results are per-assembly, most of the outputs use an assembly ID as a prefix. For the purposes of this documentation, we are using the assembly ID `ERZ12345`.

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

- **ERZ12345_filtered_contigs.fasta.gz**: This `FASTA` file contains the filtered contigs after the removal of those that are shorter than 500 bases, and which have a proportion of ambiguous bases higher than 10%.
- **ERZ12345.tsv**: This `tsv` file contains the QUAST summary output, giving an assessment of the quality of the contigs of this assembly.
- **multiqc_report.html**: This `html` file contains the `MultiQC` report for that assembly. It combines outputs from multiple tools, including `QUAST` (run both before and after quality control during assembly preprocessing), as well as records of the software versions used by the pipeline.
- **multiqc_data/**: This `directory` contains the input files used by MultiQC to generate its report.

### cds

The `cds` directory contains output files related to the combined gene caller subworkflow, which calls genes using both `Pyrodigal` and `FragGeneScanRs`, which are Python interface to Prodigal and a Rust re-implementation of FragGeneScan, respectively. The identified genes are then merged, and the resulting outputs are saved in this directory in three different formats.

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

- **ERZ12345_predicted_orf.ffn.gz**: This `FASTA` file contains the nucleotide sequences of each predicted Open Reading Frame (ORF) by the combined gene caller, making up the collection of genes detected in the assembly.
- **ERZ12345_predicted_cds.faa.gz**: This `FASTA` file contains the amino acid sequences of each predicted ORF, making up the collection of proteins detected in the assembly.
- **ERZ12345_predicted_cds.gff.gz**: This `gff` file contains the Coding DNA Sequences regions in the GFF3 format.

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
- **ERZ12345_SSU.fasta.gz**: This `FASTA` file contains all sequences from the assembly's contigs that matched the `SSU` (Small Subunit) rRNA marker gene, as identified by running `cmsearch`. These sequences represent regions of the contigs that align to the `SSU` model used in the search. This file may be absent if no marker genes of this type were detected in the assembly.
- **ERZ12345_LSU.fasta.gz**: This `FASTA` file contains sequences matching the `LSU` (Large Subunit) rRNA marker gene, identified in the same way via `cmsearch`. This file may be absent if no marker genes of this type were detected in the assembly.

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

This subdirectory contains the output generated by running `InterProScan` on the protein sequences of the assembly.

- **ERZ12345_interproscan.tsv.gz**: This `tsv` file contains the detailed results from `InterProScan`, including the InterPro identifiers, descriptions of the functional domains or protein families detected, and the specific proteins from the assembly that match each signature. Each entry links a protein to its associated functional features, such as protein families, motifs, domains, and Gene Ontology (GO) terms.
- **ERZ12345_interpro_summary.tsv.gz**: This `tsv` file provides a summary of the count of different InterPro signatures that were detected across all proteins in the assembly. The summary includes the total number of proteins annotated with each InterPro signature, offering a high-level overview of the functional diversity present in the assembly.
- **ERZ12345_interpro_summary.tsv.gz.gzi**: This file is an index for the InterPro count summary file.

#### Output files - pfam

This subdirectory contains summaries of Pfam domain matches extracted from the `InterProScan` results.

- **ERZ12345_pfam_summary.tsv.gz**: This `tsv` file contains a summary of the counts of the different Pfam signatures identified across all proteins in the assembly. Each entry represents a specific Pfam family or domain, with associated counts indicating how many proteins in the assembly match the signature.
- **ERZ12345_pfam_summary.tsv.gz.gzi**: This file is an index for the Pfam count summary file.

#### Output files - go

This subdirectory contains output files summarising Gene Ontology (GO) annotations extracted from the InterProScan results. In addition to the full set of GO terms, the annotations are also mapped to the metagenomics GO-slim subset ([see here](https://geneontology.org/docs/go-subset-guide/)) to provide a more concise and domain-specific overview.

- **ERZ12345_go_summary.tsv.gz**: This `tsv` file contains summary counts of the full set of GO terms assigned to proteins in the assembly.
- **ERZ12345_go_summary.tsv.gz.gzi**: This file is an index for the GO summary count file.
- **ERZ12345_goslim_summary.tsv.gz**: This `tsv` file contains summary counts of the GO-slim terms (a simplified subset of GO terms) assigned to proteins in the assembly.
- **ERZ12345_goslim_summary.tsv.gz.gzi**: This file is an index for the GO-slim summary count file.

#### Output files - eggnog

This subdirectory contains the output of running `EggNOG-mapper` on the protein sequences from the assembly.

- **ERZ12345_emapper_seed_orthologs.tsv.gz**: This `tsv` file lists the seed orthologs identified for each protein during the `EggNOG-mapper` run.
- **ERZ12345_emapper_annotations.tsv.gz**: This `tsv` file contains the functional annotations assigned to the proteins by `EggNOG-mapper` using the seed orthologs.

#### Output files - kegg

This subdirectory contains the results of identifying KEGG Orthologs (KOs) in the assembly’s proteins using `hmmsearch`, with their occurrences summarised in a count table.

- **ERZ12345_ko_summary.tsv.gz**: This `tsv` file listing the counts of each KO detected in the assembly.
- **ERZ12345_ko_summary.tsv.gz.gzi**: This file is an index for the KO summary count file.

#### Output files - rhea-reactions

This subdirectory contains results from running `DIAMOND` to match the assembly’s proteins against a UniRef90 database with Rhea reactions and corresponding ChEBI compounds assigned.

- **ERZ12345_proteins2rhea.tsv.gz**: This `tsv` file listing the Rhea reaction IDs and corresponding ChEBI compounds assigned to the proteins in the assembly based on their matches to UniRef90 database.
- **ERZ12345_proteins2rhea.tsv.gz.gzi**: This file is an index for the Proteins2Rhea output.

#### Output files - dbcan

This subdirectory contains the output of running `run_dbCAN` on the assembly proteins, generating annotations for CAZymes and CAZyme Gene Clusters (CGC).

- **ERZ12345_dbcan_cgc.gff.gz**: This `gff` file is annotated with functional genes for CGCFinder and visualization tools.
- **ERZ12345_dbcan_standard_out.tsv.gz**: This `tsv` file lists all identified CGCs and their components.
- **ERZ12345_dbcan_overview.tsv.gz**: This `tsv` file contains a summary of identified CAZymes.
- **ERZ12345_dbcan_hmm.tsv.gz**: This `tsv` file contains the detailed HMMER results.
- **ERZ12345_dbcan_sub_hmm.tsv.gz**: This `tsv` file contains the detailed sub-HMMER results.
- **ERZ12345_dbcan_substrates.tsv.gz**: This `tsv` file contains the substrate prediction results for CGCs.

### pathways-and-systems

The `pathways-and-systems` directory contains five subdirectories of results, one for each pathway/system annotation tool used by the pipeline. The results range from KEGG pathways, to Biosynthetic Gene Clusters (BGCs) detected by `antiSMASH` and `SanntiS`. Just like the functional annotations, most of these annotations are on a per-protein basis.

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

This subdirectory contains the results of running `AntiSMASH` on the assembly’s proteins, including annotations for detected BGCs in various formats.

- **ERZ12345_antismash.gbk.gz**: A GenBank format file containing AntiSMASH annotations.
- **ERZ12345_antismash.gff.gz**: A GFF3 format file containing AntiSMASH annotations.
- **ERZ12345_merged.json**: A JSON file containing AntiSMASH annotations.
- **ERZ12345_antismash_summary.tsv.gz**: This `tsv` file listing the counts of each AntiSMASH entry detected in the assembly.

#### Output files - sanntis

This subdirectory contains the outputs of running `SanntiS` on the proteins of the assembly, also describing the detected BGCs using this new machine learning-based tool.

- **ERZ12345_sanntis.gff.gz**: This `gff` file contains the different SanntiS annotations in the GFF3 format. The GFF3 specification can be found in the [tool repo](https://github.com/Finn-Lab/SanntiS?tab=readme-ov-file#ouput)

#### Output files - genome-properties

This subdirectory contains the outputs of running `Genome Properties` on the proteins of the assembly. Genome properties is an annotation system whereby functional attributes can be assigned to a assembly, based on the presence of a defined set of protein signatures within that assembly. The results are available in three different formats.

- **ERZ12345_gp.json.gz**: A JSON format file containing Genome Properties annotations
- **ERZ12345_gp.tsv.gz**: A `tsv` file containing Genome Properties annotations in a tab-separated format.
- **ERZ12345_gp.txt.gz**: A plain text file containing the Genome Properties annotations.

#### Output files - kegg-modules

This subdirectory contains results from evaluating which KEGG modules and pathways are represented in the assembly based on the presence of KEGG Orthologs (KOs). For each module, a completeness score is calculated depending on how many required KOs are detected using the [kegg-pathways-completeness-tool](https://github.com/EBI-Metagenomics/kegg-pathways-completeness-tool).

- **ERZ12345_kegg_modules_per_contigs.tsv.gz**: This `tsv` file contains the different KEGG modules and their completeness, including the contig they were found on.
- **ERZ12345_kegg_modules_per_contigs.tsv.gz.gzi**: This file is an index for the per-contig KEGG modules file.
- **ERZ12345_kegg_modules_summary.tsv.gz**: This `tsv` file contains a summary of the KEGG modules and their completeness for the whole assembly, not on a per-contig basis.
- **ERZ12345_kegg_modules_summary.tsv.gz.gzi**: This file is an index for the summary KEGG modules file.

#### Output files - dram-distill

This subdirectory contains the outputs of running `DRAM-distill` on the KO hits per contig, the InterProScan and run_dbCan. DRAM generates summary reports and visualisations for the annotations on a per-assembly basis, and outputs four different files.

- **ERZ12345_dram.tsv.gz**: This `tsv` file contains the product of `DRAM-distill` for the assembly, including functions that were detected.
- **ERZ12345_dram.html.gz**: This `html` file contains a heatmap visualisation of the detected functions by `DRAM-distill`.

### gff

The `gff` directory contains a single file that summarises the annotations generated by the pipeline, integrating functional annotation per protein and other genomic features.

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

- **ERZ12345_annotation_summary.gff.gz**: This `gff` file contains an expansive and large summary of most of the functional annotations each protein has into a single GFF3 format file.

## Per-study output files

The pipeline generates four different per-study output files that aggregate and summarise data, from successful assembly analysis runs to study-wide `MultiQC` reports. These are stored at the root of the study analysis directory.

### Successfully analysed assemblies

The IDs of assemblies that have been successfully analysed are aggregated into a top-level file (`analysed_assemblies.csv`), which looks like this:

```
ERZ12345,success
ERZ56789,success
```

### MultiQC

In addition to generating MultiQC reports for each assembly, the pipeline also produces per-study reports, which can be found in the `multiqc/` directory.

### dram-distill

Just as the pipeline generates `DRAM-distill` outputs on a per-assembly basis, it also produces similar outputs on a per-study basis, which are located in the `dram-distill/` directory by concatenating assembly-level files.
