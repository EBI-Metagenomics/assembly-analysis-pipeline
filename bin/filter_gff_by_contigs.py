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
import fileinput
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def contig_ids_from_fasta(fasta_file: Path) -> set[str]:
    """Read the sequence IDs of a FASTA file.

    Only the first word of each header is used, the same as seqkit and antiSMASH.

    :param fasta_file: The FASTA file. It can be compressed (``.gz``).
    :return: The contig IDs.
    """
    contig_ids: set[str] = set()
    with fileinput.input(fasta_file, openhook=fileinput.hook_compressed) as fasta:
        for line in fasta:
            if line.startswith(">"):
                contig_ids.add(line[1:].split()[0])
    return contig_ids


def filter_gff_by_contigs(
    gff_file: Path, contig_ids: set[str], output_file: Path
) -> int:
    """Write the GFF directives (``##``) and the features of the given contigs.

    Comment lines (``#``) are not written, because the CGC GFF has one comment line
    for each contig.

    :param gff_file: The GFF file. It can be compressed (``.gz``).
    :param contig_ids: The contigs to keep (column 1 of the GFF).
    :param output_file: The output GFF file. It is compressed if the name ends in ``.gz``.
    :return: The number of CDS features in the output.
    """
    cds_count = 0
    with (
        fileinput.input(gff_file, openhook=fileinput.hook_compressed) as gff,
        fileinput.hook_compressed(output_file, "w") as output,
    ):
        for line in gff:
            if line.startswith("##"):
                output.write(line)
                continue
            if line.startswith("#") or not line.strip():
                continue

            columns = line.split("\t")
            if columns[0] not in contig_ids:
                continue

            output.write(line)
            if len(columns) > 2 and columns[2] == "CDS":
                cds_count += 1
    return cds_count


def main() -> None:
    """Filter a GFF to the contigs of a FASTA file.

    antiSMASH reads the full GFF for each contigs chunk. With this filter, each chunk
    gets only its own features. antiSMASH fails if no CDS is in the GFF, so the output
    file is removed if it has no CDS.
    """
    parser = argparse.ArgumentParser(
        description="Filter a GFF to the contigs of a FASTA file."
    )
    parser.add_argument(
        "--contigs",
        required=True,
        type=Path,
        help="Contigs FASTA file (optionally gzipped).",
    )
    parser.add_argument(
        "--gff",
        required=True,
        type=Path,
        help="GFF file to filter (optionally gzipped).",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output GFF file (gzipped if the name ends in .gz).",
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    contig_ids = contig_ids_from_fasta(args.contigs)
    logger.info(f"Found {len(contig_ids)} contigs in '{args.contigs}'.")

    cds_count = filter_gff_by_contigs(args.gff, contig_ids, args.output)
    if cds_count == 0:
        logger.warning(
            f"No CDS for the contigs in '{args.contigs}'. No output is written."
        )
        args.output.unlink()
        return

    logger.info(f"Wrote {cds_count} CDS to '{args.output}'.")


if __name__ == "__main__":
    main()
