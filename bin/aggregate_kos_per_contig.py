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

import csv
from collections import defaultdict
import argparse
import gzip


def convert_ko_per_contig(input_file, output_file):
    """
    Converts the input file with 'ko' and 'contig_id' columns
    into a tab-separated file where each row contains a contig ID
    followed by all the associated KO IDs for that contig.

    For example, the input file:
    ko      contig_id
    K07082  MGYA1767_1
    K09458  MGYA1767_1
    K02990  MGYA1767_1

    will be converted to the output file:
    MGYA1767_1    K07082  K09458  K02990
    """
    ko_per_contig = defaultdict(set)

    with gzip.open(input_file, "rt") as file_handle:
        reader = csv.DictReader(
            file_handle, delimiter="\t", fieldnames=["ko", "contig_id"]
        )
        for row in reader:
            ko_per_contig[row["contig_id"]].add(row["ko"])

    with open(output_file, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        for contig_id, kos in ko_per_contig.items():
            writer.writerow([contig_id] + list(kos))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate KO IDs per contig")
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help='Input file path (tab-separated with "ko" and "contig_id" columns)',
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output file path (tab-separated with contig ID followed by associated KO IDs)",
    )
    args = parser.parse_args()

    convert_ko_per_contig(args.input, args.output)
