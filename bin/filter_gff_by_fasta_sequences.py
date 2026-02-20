#!/usr/bin/env python3

import logging
import os
import re
import sys
import argparse


logger = logging.getLogger(__name__)

# Prodigal ID token: SAMPLE_CONTIG-NUM_GENE-NUM (last two parts are integers)
_PRODIGAL_ID_RE = re.compile(r"^.+_\d+_\d+$")


def _parse_fasta_header(header: str, lineno: int) -> str:
    """Parse a single FASTA header line and return the corresponding sequence ID.

    Detects whether the header belongs to FragGeneScanRS or Prodigal and
    applies format-specific validation before returning the ID.

    :param header: Raw header line, including the leading ``>``.
    :param lineno: Line number in the source file, used in error messages.
    :returns: Extracted sequence identifier.
    """
    token = header[1:].split()[0]
    parts = token.split("_")

    # FragGeneScanRS detection: last field must be a strand indicator #
    if parts[-1] in {"+", "-"}:
        if len(parts) < 4:
            raise ValueError(
                f"Line {lineno}: FragGeneScanRS header '{token}' has fewer than "
                f"4 underscore-separated fields (got {len(parts)})."
            )
        return token

    # Prodigal detection: ID token must match SAMPLE_CONTIG_GENE #
    if _PRODIGAL_ID_RE.match(token):
        return token

    raise ValueError(
        f"Line {lineno}: unrecognised header format '{token}'. "
        f"Expected Prodigal (e.g. 'ERZ23524530_111_2') or "
        f"FragGeneScanRS (e.g. 'ERZ23524530_383157_129_194_-')."
    )


def _parse_gff_id(attributes: str, lineno: int, gff_file: str) -> str:
    """Extract the ``ID=`` value from a GFF attributes field (column 9).

    :param attributes: Raw attributes string (e.g. ``ID=ERZ102_1_1;other=val``).
    :param lineno: Line number in the source file, used in error messages.
    :param gff_file: Path to the GFF file, used in error messages.
    :returns: The value of the ``ID`` attribute.
    :raises ValueError: If no ``ID=`` field is found in the attributes.
    """
    for field in attributes.split(";"):
        if field.startswith("ID="):
            return field[3:]
    raise ValueError(
        f"Malformed GFF line {lineno} in '{gff_file}': "
        f"no ID= field found in attributes '{attributes}'."
    )


def extract_fasta_sequence_ids(fasta_file: str) -> set[str]:
    """Extract unique sequence IDs from a FASTA file.

    :param fasta_file: Path to the input FASTA file.
    :returns: Set of extracted sequence identifiers.
    """
    fasta_ids: set[str] = set()

    with open(fasta_file, "r") as fh:
        for lineno, line in enumerate(fh, start=1):
            if not line.startswith(">"):
                continue
            seq_id = _parse_fasta_header(line.rstrip("\n"), lineno)
            fasta_ids.add(seq_id)

    return fasta_ids


def filter_gff_by_sequence_ids(
    gff_file: str, sequence_ids: set[str], output_file: str
) -> set[str]:
    """Write a filtered GFF retaining only features matching *sequence_ids*.

    Matching is done against the ``ID=`` value in the attributes column (column 9).
    Header lines (starting with ``#``) are always preserved.

    :param gff_file: Path to the input GFF file.
    :param sequence_ids: Set of sequence IDs to retain.
    :param output_file: Path for the filtered output GFF file.
    :returns: Set of all gene IDs encountered in *gff_file*.
    """
    gff_sequence_ids: set[str] = set()

    with open(gff_file, "r") as infile, open(output_file, "w") as outfile:
        for line_num, line in enumerate(infile, start=1):

            if line.startswith("#"):
                outfile.write(line)
                continue

            if not line.strip():
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                raise ValueError(
                    f"Malformed GFF line {line_num} in '{gff_file}': "
                    f"expected >=9 fields, got {len(parts)}."
                )
            seq_id = _parse_gff_id(parts[8], line_num, gff_file)
            gff_sequence_ids.add(seq_id)
            if seq_id in sequence_ids:
                outfile.write(line)

    return gff_sequence_ids


def main() -> None:
    """Filter a GFF file to retain only features whose sequence IDs are present in a FASTA file.

    Supports two gene-calling output formats:

    - Prodigal: ``>ERZ23524530_111_2 # 841 # 11490 # -1 # ID=12368_2;...``
    - FragGeneScanRS: ``>ERZ23524530_383157_129_194_-``

    For FragGeneScanRS entries the base contig ID is reconstructed by keeping only
    the first two underscore-delimited fields (e.g. ``ERZ23524530_383157``).
    """
    parser = argparse.ArgumentParser(
        description="Filter a GFF file to only include features for sequences present in a FASTA file."
    )
    parser.add_argument(
        "fasta_file", help="Input FASTA file containing protein sequences."
    )
    parser.add_argument("gff_file", help="Input GFF file to be filtered.")
    parser.add_argument("output_file", help="Output filtered GFF file.")
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose (DEBUG-level) logging.",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    if not os.path.isfile(args.fasta_file):
        logger.error(f"FASTA file '{args.fasta_file}' not found.")
        sys.exit(1)

    if not os.path.isfile(args.gff_file):
        logger.error(f"GFF file '{args.gff_file}' not found.")
        sys.exit(1)

    sequence_ids = extract_fasta_sequence_ids(args.fasta_file)

    logger.info(
        f"Found {len(sequence_ids)} unique sequence IDs in '{args.fasta_file}'."
    )

    gff_sequence_ids = filter_gff_by_sequence_ids(
        args.gff_file, sequence_ids, args.output_file
    )

    missing = sequence_ids - gff_sequence_ids
    if missing:
        missing_list = ", ".join(sorted(missing))
        logger.error(
            f"{len(missing)} FASTA sequence(s) are absent from the GFF file: {missing_list}"
        )
        sys.exit(1)

    logger.info("Done.")


if __name__ == "__main__":
    main()
