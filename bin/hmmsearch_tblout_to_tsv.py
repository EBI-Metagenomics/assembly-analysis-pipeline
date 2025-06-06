#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2025 EMBL - European Bioinformatics Institute
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

import sys
import csv
from Bio import SearchIO
import argparse
import gzip


def parse_hmmsearch(input_file, output_handle):
    with gzip.open(input_file, "rt") as input_handle:
        writer = csv.writer(output_handle, delimiter="\t")
        for qresult in SearchIO.parse(input_handle, "hmmer3-tab"):
            ko_id = qresult.id
            for hit in qresult.hits:
                writer.writerow([ko_id, hit.id])


if __name__ == "__main__":
    """
    This script parses the output of the hmmsearch - tblout and writes a subset as a TSV

    It reads from the input file specified as an argument and writes to standard output.

    The output format is: KO_ID<tab>CONTIG_ID
    """
    parser = argparse.ArgumentParser(
        description="Parse hmmsearch tblout and write a subset as TSV"
    )
    parser.add_argument("input_file", help="Input file path (gzipped)")
    args = parser.parse_args()

    parse_hmmsearch(args.input_file, sys.stdout)
