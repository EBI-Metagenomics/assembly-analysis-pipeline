# Agent initialization — assembly-analysis-pipeline

This file provides context for AI coding agents working in this repository.

## Project overview

nf-core DSL2 Nextflow pipeline for metagenome assembly annotation (MGnify v6+).
Repository: `ebi-metagenomics/assembly-analysis-pipeline`

## Coding conventions

- Nextflow DSL2 only, following nf-core conventions.
- Unit tests use [nf-test](https://www.nf-test.com/). Every module/subworkflow change should have a corresponding test.
- Python code uses type hints, reST-style docstrings, and is formatted with black.
- Do not run Python scripts directly — the environment is Dockerized. All tools run inside containers.
- Use `gh` (GitHub CLI) for any GitHub interactions.

## Key file paths

| Path                                      | Purpose                                       |
| ----------------------------------------- | --------------------------------------------- |
| `main.nf`                                 | Entrypoint                                    |
| `workflows/assembly_analysis_pipeline.nf` | Main workflow logic                           |
| `subworkflows/local/`                     | Pipeline-specific subworkflows                |
| `subworkflows/ebi-metagenomics/`          | Shared EBI-Metagenomics subworkflows          |
| `modules/local/`                          | Pipeline-local modules                        |
| `modules/ebi-metagenomics/`               | Shared EBI-Metagenomics modules               |
| `modules/nf-core/`                        | Imported nf-core modules                      |
| `conf/base.config`                        | Resource label definitions                    |
| `conf/modules.config`                     | Per-process publishDir and argument overrides |
| `nextflow_schema.json`                    | Parameter schema                              |
| `assets/samplesheet.csv`                  | Example samplesheet                           |

## Pipeline stages (in order)

1. `DOWNLOAD_FROM_FIRE` — optional EBI-internal FIRE S3 download
2. `RENAME_CONTIGS` — standardise contig names to `{id}_{n}`
3. `ASSEMBLY_QC` — length/N-content filtering, optional decontamination (minimap2), QUAST
4. `DETECT_RNA` — Infernal cmsearch, deoverlap, extract SSU/LSU sequences
5. `COMBINED_GENE_CALLER` — Pyrodigal + FragGeneScanRS, merged into GFF/FAA/FFN
6. `CONTIGS_TAXONOMIC_CLASSIFICATION` — CAT_pack + Krona
7. `SEQKIT_SPLIT2` — chunk proteins for parallel annotation
8. `FUNCTIONAL_ANNOTATION` — InterProScan, eggNOG-mapper, dbCAN, KEGG KOs (hmmsearch), RHEA (Diamond), GO/GOslim, Pfam summaries
9. `PATHWAYS_AND_SYSTEMS` — SANNTIS, antiSMASH, KEGG modules, Genome Properties, DRAM-distill
10. `GFF_SUMMARY` — consolidated annotation GFF3
11. `GT_GFF3VALIDATOR` — GFF validation (non-fatal; errors written to log)
12. `MULTIQC` — per-assembly and per-samplesheet reports

## Resource labels (conf/base.config)

All values scale by `task.attempt`. Default `maxRetries` is 1.

| Label                 | CPUs | Memory | Time |
| --------------------- | ---- | ------ | ---- |
| `process_single`      | 1    | 6 GB   | 4 h  |
| `process_low`         | 2    | 12 GB  | 4 h  |
| `process_medium`      | 6    | 36 GB  | 8 h  |
| `process_high`        | 12   | 72 GB  | 16 h |
| `process_high_memory` | —    | 200 GB | —    |
| `process_long`        | —    | —      | 20 h |

Resource overrides for specific processes go in `conf/modules.config` using `withName`.

## Output structure

Each assembly (e.g. `ERZ12345`) produces outputs under `{outdir}/{id}/`:

```
qc/
cds/
taxonomy/
functional-annotation/
  interpro/  pfam/  go/  eggnog/  kegg/  rhea-reactions/  dbcan/
pathways-and-systems/
  antismash/  sanntis/  genome-properties/  kegg-modules/  dram-distill/
annotation-summary/
```

Study-level outputs at `{outdir}/`: `analysed_assemblies.csv`, `qc_failed_assemblies.csv`,
`multiqc/`, `dram-distill/`, `downstream_samplesheets/`, `pipeline_info/`.

## Samplesheet format

```csv
sample,assembly_fasta,contaminant_reference,human_reference,phix_reference
ERZ999,/path/to/ERZ999.fasta.gz,,,
```

The decontamination columns are optional. References are subfolder names within
`--reference_genomes_folder`, each containing a `.fna` file.
