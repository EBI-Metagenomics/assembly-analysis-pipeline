# Python unit tests

Unit tests for scripts in `bin/`. Uses the standard `unittest` module — no extra dependencies required.

## Run all tests

```bash
python -m unittest discover -s tests/python -v
```

## Run a single test file

```bash
python tests/python/filter_ips_by_contigs_test.py -v
python tests/python/filter_gff_by_fasta_sequences_test.py -v
```

## Run a single test case or method

```bash
python tests/python/filter_ips_by_contigs_test.py TestContigIdFromProteinId -v
python tests/python/filter_ips_by_contigs_test.py TestContigIdFromProteinId.test_fgs_format_negative_strand -v
```

## Test coverage

| Test file                               | Script under test                      | Coverage                                                                                                                                |
| --------------------------------------- | -------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------- |
| `filter_ips_by_contigs_test.py`         | `bin/filter_ips_by_contigs.py`         | FASTA parsing, protein ID resolution (Pyrodigal and FragGeneScanRS formats), IPS TSV filtering, gzip output, stdout protein ID emission |
| `filter_gff_by_fasta_sequences_test.py` | `bin/filter_gff_by_fasta_sequences.py` | FASTA sequence ID extraction, GFF filtering by sequence ID, header preservation                                                         |
