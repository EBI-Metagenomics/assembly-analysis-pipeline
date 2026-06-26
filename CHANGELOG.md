# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [6.0.3] - 2026-06-26

### Changed

- Increased the runtime and memory minimap2_align uses by default.

## [6.0.2] - 2026-04-24

### Fixed

- Fixed RHEA annotation output to parse `contig_id` values correctly from protein identifiers and standardised the output header and `top_hit` values.

## [6.0.1] - 2026-04-15

### Added

- Added repository instructions for coding agents in `AGENTS.md`.
- Added an `nf-test` covering EggNOG annotation TSV concatenation.

### Changed

- Switched merged EggNOG annotation TSV generation to `csvtk concat` so the combined output keeps a single preserved header and remains compressed as `tsv.gz`.

### Fixed

- Fixed duplicated header rows appearing in concatenated EggNOG annotation TSV outputs.

## [6.0.0] - 2026-03-13

### Added

- First stable release of the MGnify metagenome assembly analysis pipeline.
- Added support for optional FIRE downloads, contig renaming, assembly QC, and optional decontamination before downstream analysis.
- Added RNA detection, combined gene calling, contig taxonomic classification, and per-assembly taxonomy visualisation.
- Added the functional annotation workflow covering InterProScan, EggNOG-mapper, dbCAN, KEGG orthologs, RHEA reactions, GO/GOslim, and Pfam summaries.
- Added pathways and systems analyses with SANNTIS, antiSMASH, KEGG modules, Genome Properties, and DRAM distillation.
- Added consolidated annotation summaries, GFF3 validation, MultiQC reporting, and downstream samplesheet generation.

### Changed

- Promoted the pipeline from the `6.0.0` beta and release-candidate series to the first stable `6.0.0` release.
- Reworked chunking and resource tuning across functional annotation and pathway modules to improve runtime and memory behaviour on large assemblies.
- Updated core tooling and shared modules throughout the release-candidate cycle, including InterProScan, antiSMASH 8, dbCAN/easySubstrate, and the combined gene caller.

### Fixed

- Fixed multiple edge cases in annotation summarisation and GFF generation, including empty-output handling and no-hit cases.
- Fixed downstream workflows to consistently use QC-filtered and decontaminated contigs.
- Fixed reporting and published outputs across MultiQC, downstream samplesheets, and compressed summary files.

[unreleased]: https://github.com/ebi-metagenomics/assembly-analysis-pipeline/compare/v6.0.2...HEAD
[6.0.2]: https://github.com/ebi-metagenomics/assembly-analysis-pipeline/compare/v6.0.1...HEAD
[6.0.1]: https://github.com/ebi-metagenomics/assembly-analysis-pipeline/compare/v6.0.0...060b74274e5a3fcf92fdf76040f95770a5f848c6
[6.0.0]: https://github.com/ebi-metagenomics/assembly-analysis-pipeline/releases/tag/v6.0.0
