#!/usr/bin/env python
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
import csv
import gzip
import sys
from typing import TextIO


def parse_hmmsearch(input_file: str, output_handle: TextIO) -> None:
    """Stream reported domain hits as KO/protein pairs.

    Comments and blank lines may occur anywhere, including between complete
    tables in concatenated gzip members. A table's footer does not end parsing.

    Example header and hit (first six whitespace-separated columns)::

        # target name  accession  tlen  query name  accession  qlen ...
          protein_1    -          300   K00001      -          100  ...

    Column 1 is the target/protein ID; column 4 is the query/KO ID.
    All reported domains are emitted. The summary step deduplicates pairs with
    csvtk before counting proteins per KO, matching v5 KO summaries.
    The example hit produces ``K00001<tab>protein_1`` in the output.

    :param input_file: Path to the gzipped HMMER ``--domtblout`` tables.
    :param output_handle: Text stream for the headerless two-column TSV.
    """
    with gzip.open(input_file, "rt") as input_handle:
        writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
        for line_number, line in enumerate(input_handle, start=1):
            # Only columns 1 (protein) and 4 (KO) are needed.
            fields = line.split(maxsplit=4)
            if not fields or fields[0].startswith("#"):
                continue
            if len(fields) < 4:
                raise ValueError(f"Invalid HMMER domtblout row at line {line_number}")
            writer.writerow([fields[3], fields[0]])


if __name__ == "__main__":
    """
    This script parses hmmsearch --domtblout and streams KO/protein pairs as TSV.

    It reads from the input file specified as an argument and writes to standard output.

    The output format is: KO_ID<tab>PROTEIN_ID
    """
    parser = argparse.ArgumentParser(
        description="Parse hmmsearch domtblout and stream KO/protein pairs as TSV"
    )
    parser.add_argument("input_file", help="Input file path (gzipped)")
    args = parser.parse_args()

    parse_hmmsearch(args.input_file, sys.stdout)
