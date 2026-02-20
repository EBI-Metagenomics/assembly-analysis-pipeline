# Tests for `filter_gff_by_fasta_sequences.py`

## Run

```bash
# All tests
python test_filter_gff.py -v

# Single test
python test_filter_gff.py TestExtractFastaSequenceIds.test_mixed_formats -v
```

## Coverage

- `TestExtractFastaSequenceIds` — Prodigal and FragGeneScanRS header parsing, including invalid formats
- `TestFilterGffBySequenceIds` — filtering logic, header preservation, malformed GFF detection
