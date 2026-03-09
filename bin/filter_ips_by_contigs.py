#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright 2026 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import argparse
import gzip
import logging
import sys
from pathlib import Path
from typing import Iterator

logging.basicConfig(
    format="%(levelname)s: %(message)s",
    level=logging.INFO,
    stream=sys.stderr,
)


def _open(path: Path | str):
    """Open a plain or gzip-compressed text file for reading.

    :param path: Path to the file.
    :return: A text-mode file object.
    """
    if Path(path).suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "r")


def _contig_ids_from_fasta(fasta_path: Path) -> set[str]:
    """Parse a FASTA file and return the set of sequence IDs.

    Only the first whitespace-delimited token of each header line is used,
    matching the behaviour of seqkit and other common tools.

    :param fasta_path: Path to the FASTA file (plain or gzip).
    :return: Set of contig IDs (without the leading '>').
    """
    ids: set[str] = set()
    with _open(fasta_path) as handle:
        for line in handle:
            if line.startswith(">"):
                ids.add(line[1:].split()[0])

    if not ids:
        logging.error(f"No sequences found in {fasta_path}")
        sys.exit(1)

    return ids


def _contig_id_from_protein_id(protein_id: str, contig_ids: set[str]) -> str | None:
    """Derive the contig ID from a protein ID, handling both CGC callers.

    The CGC produces two protein ID formats:

    - Pyrodigal: ``{contig_id}_{cds_n}`` — strip the last component.
    - FragGeneScanRS: ``{contig_id}_{start}_{end}_{strand}`` — strip the last three.

    Both are tried in order; the first match against *contig_ids* is returned.
    Returns ``None`` if neither resolves to a known contig.

    :param protein_id: Protein identifier from the IPS TSV.
    :param contig_ids: Set of known contig IDs.
    :return: Matching contig ID, or ``None``.
    """
    for candidate in (protein_id.rsplit("_", 1)[0], protein_id.rsplit("_", 3)[0]):
        if candidate in contig_ids:
            return candidate
    return None


def _iter_filtered_rows(
    ips_path: Path,
    contig_ids: set[str],
) -> Iterator[tuple[str, str]]:
    """Stream IPS TSV rows that belong to contigs in *contig_ids*.

    The protein ID is column 1 (0-indexed: column 0) of the tab-separated file.
    Both Pyrodigal (``{contig_id}_{cds_n}``) and FragGeneScanRS
    (``{contig_id}_{start}_{end}_{strand}``) protein ID formats are handled.

    :param ips_path: Path to the InterProScan TSV (plain or gzip).
    :param contig_ids: Set of contig IDs to keep.
    :return: Iterator of ``(protein_id, raw_line)`` tuples for matching rows.
    :raises SystemExit: If a malformed row (no tab separator) is encountered.
    """
    with _open(ips_path) as handle:
        for lineno, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t", 1)
            if len(parts) < 2:
                logging.error(
                    f"Malformed IPS row at line {lineno} (no tab separator): {line.rstrip()}"
                )
                sys.exit(1)
            protein_id = parts[0]
            contig_id = _contig_id_from_protein_id(protein_id, contig_ids)
            if contig_id is not None:
                yield protein_id, line


def filter_ips(
    contigs_path: Path,
    ips_path: Path,
    out_path: Path,
) -> None:
    """Filter *ips_path* to rows whose proteins belong to the provided contigs.

    Each matched protein ID is written to stdout (one per line) so the caller
    can pipe them into ``seqkit grep -f -``.

    :param contigs_path: FASTA file containing the contig subset.
    :param ips_path: Full InterProScan TSV to filter.
    :param out_path: Destination for the filtered TSV (written as gzip).
    """
    contig_ids = _contig_ids_from_fasta(contigs_path)

    with gzip.open(out_path, "wt") as out_fh:
        for protein_id, line in _iter_filtered_rows(ips_path, contig_ids):
            out_fh.write(line)
            sys.stdout.write(protein_id + "\n")
    sys.stdout.flush()


if __name__ == "__main__":
    """Filter an InterProScan TSV to rows whose proteins belong to a given set of contigs.

    Matched protein IDs are written to stdout so the caller can pipe them directly
    into ``seqkit grep -f -`` to subset a FAA file without creating a temporary file.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--contigs",
        required=True,
        type=Path,
        metavar="FASTA",
        help="Contig chunk FASTA (plain or gzip).",
    )
    parser.add_argument(
        "--ips",
        required=True,
        type=Path,
        metavar="TSV",
        help="Full InterProScan TSV (plain or gzip).",
    )
    parser.add_argument(
        "--out",
        required=True,
        type=Path,
        metavar="TSV_GZ",
        help="Output filtered TSV (written as gzip).",
    )
    args = parser.parse_args()
    filter_ips(
        contigs_path=args.contigs,
        ips_path=args.ips,
        out_path=args.out,
    )
