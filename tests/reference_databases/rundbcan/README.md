# Tiny dbCAN database (test fixture)

A minimal subset of the [run_dbcan](https://github.com/bcb-unl/run_dbcan) database, used by
`conf/test.config` / `conf/test_full.config` (`--dbcan_database`). It keeps the exact
file/folder layout and file formats of the full database, but with only a handful of entries,
so `run_dbcan easy_substrate` runs end-to-end in seconds.

Built from the full database used by `dbcan=5.2.9` (`quay.io/biocontainers/dbcan:5.2.9--pyhdfd78af_0`).

## Contents

All entries are centred on **PUL0111** (_E. coli_ K-12 melibiose locus: melR / melA / melB),
which matches the test proteins in `ecoliK12MG1655_test.faa`.

| File                            | Entries kept                                                           |
| ------------------------------- | ---------------------------------------------------------------------- |
| `dbCAN.hmm`                     | `GH4.hmm`                                                              |
| `dbCAN-sub.hmm`                 | `GH4_e0.hmm\|GH4:2933\|3.2.1.22:13`                                    |
| `STP.hmm`                       | `Fer4`, `HTH_18`                                                       |
| `TF.hmm`                        | `HLH`, `HTH_3`                                                         |
| `CAZy.dmnd`                     | `AAC77080.1` (melA), `QSR40788.1`                                      |
| `PUL.dmnd`                      | `PUL0111_1..3`                                                         |
| `TCDB.dmnd`                     | 35 TC-DB entries (incl. `A7ZUZ0`, melB) — list in the script           |
| `TF.dmnd`                       | `P0ACH8` (MelR), `P0A9E0` (AraC)                                       |
| `peptidase_db.dmnd`             | `MER0000001`                                                           |
| `sulfatlas_db.dmnd`             | `A5AB00_S1_4`                                                          |
| `fam-substrate-mapping.tsv`     | full copy (small)                                                      |
| `dbCAN-PUL.xlsx`                | full copy (small)                                                      |
| `dbCAN-PUL/PUL0111.out/`        | copy of the single PUL folder                                          |
| `ecoliK12MG1655_test.{faa,gff}` | test inputs (3 proteins from the mel operon region), not part of dbCAN |

## How to regenerate

Requirements: `diamond` (>= 2.1, database format v3), `awk`, `bash`, and the full dbCAN
database downloaded locally (e.g. `run_dbcan database --db_dir DBCAN`).

1. Check that the entries above still exist in the new release, as names/headers change
   between versions (e.g. `grep '^NAME' DBCAN/TF.hmm`, `diamond getseq -d DBCAN/PUL.dmnd | grep '>PUL0111_'`).
   Adjust the names in the script if needed. Note that the `dbCAN-sub.hmm` model name embeds
   sequence counts (`GH4:2933`), so it will likely change.
2. Run the script below from the repository root, with `D` pointing to the full database:

```bash
D=DBCAN                                   # full database
O=tests/reference_databases/rundbcan_new  # output
W=$(mktemp -d)
mkdir -p $O/dbCAN-PUL

# HMM files: extract whole models by exact NAME (comma-separated list)
hmmx() {
  LC_ALL=C awk -v names="$3" '
    BEGIN { n = split(names, a, ","); for (i = 1; i <= n; i++) want[a[i]] = 1 }
    { buf = buf $0 "\n" }
    /^NAME/ { nm = $0; sub(/^NAME +/, "", nm); keep = (nm in want) }
    /^\/\// { if (keep) printf "%s", buf; buf = ""; keep = 0 }' "$1" > "$2"
}
hmmx $D/dbCAN.hmm     $O/dbCAN.hmm     "GH4.hmm"
hmmx $D/dbCAN-sub.hmm $O/dbCAN-sub.hmm "GH4_e0.hmm|GH4:2933|3.2.1.22:13"
hmmx $D/STP.hmm       $O/STP.hmm       "Fer4,HTH_18"
hmmx $D/TF.hmm        $O/TF.hmm        "HLH,HTH_3"

# DIAMOND files: dump sequences, keep headers matching a regex, rebuild.
# Use [|] and [.] rather than \| and \. (BSD/macOS awk rejects those escapes).
dmx() {
  local name=$(basename "$1" .dmnd)
  diamond getseq -d "$1" 2>/dev/null \
    | LC_ALL=C awk -v pat="$2" '/^>/ { p = ($0 ~ pat) } p' > $W/$name.faa
  diamond makedb --quiet --in $W/$name.faa -d $O/$name
}
TC='A7ZUZ0|B0SM05|B1LPP9|P07658|P0AAE8|P0ABN5|P0ABN9|P0ADB7|P0AF52|P0AF54|P0C0L7|P16677|P16682|P17596|P21345|P31076|P31077|P32703|P32705|P32714|P32715|P32720|P32721|P36655|P39265|P39276|P39277|P39282|P39285|P60061|P69937|P75957|Q97VF5|Q9I3N2|X5H3P0'
dmx $D/CAZy.dmnd         '^>(AAC77080[.]1|QSR40788[.]1)[|]'
dmx $D/PUL.dmnd          '^>PUL0111_'
dmx $D/TCDB.dmnd         "^>gnl[|]TC-DB[|]($TC)[|]"
dmx $D/TF.dmnd           '^>sp[|](P0ACH8|P0A9E0)[|]'
dmx $D/peptidase_db.dmnd '^>MER0000001[|]'
dmx $D/sulfatlas_db.dmnd 'A5AB00_S1_4'

# Small files: copy as-is
cp $D/fam-substrate-mapping.tsv $D/dbCAN-PUL.xlsx $O/
cp -R $D/dbCAN-PUL/PUL0111.out $O/dbCAN-PUL/

# Keep the test inputs and this README
cp tests/reference_databases/rundbcan/{ecoliK12MG1655_test.faa,ecoliK12MG1655_test.gff,README.md} $O/
```

3. Sanity check: every HMM file has the expected models, and every `.dmnd` has the expected
   number of sequences (`diamond dbinfo -d $O/<file>.dmnd`). If the new release adds a
   new database file, add it too (in `v5` `TF.dmnd` was added).
4. Smoke test with the same container version as `modules/nf-core/rundbcan/easysubstrate`:

```bash
docker run --rm --platform linux/amd64 -v $PWD/$O:/db:ro quay.io/biocontainers/dbcan:5.2.9--pyhdfd78af_0 bash -c '
  cd /tmp && run_dbcan easy_substrate --mode protein \
    --input_raw_data /db/ecoliK12MG1655_test.faa --input_gff /db/ecoliK12MG1655_test.gff \
    --gff_type prodigal --db_dir /db --output_dir res --threads 2 \
  && cat res/overview.tsv res/cgc_standard_out.tsv res/substrate_prediction.tsv'
```

Expected: `NC_000913.3_4037` annotated as `GH4_e0` (raffinose), one CGC (CAZyme + TC),
and `substrate_prediction.tsv` matching `PUL0111` / `melibiose`.

5. Replace `tests/reference_databases/rundbcan/` with `$O` and update the nf-test snapshots
   (`nf-test test ... --update-snapshot`).
