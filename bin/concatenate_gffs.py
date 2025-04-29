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
import argparse
import gzip
from typing import List


def concatenate_gff_files(gffs: List[str], output_file: str) -> None:
    with gzip.open(output_file, "wt", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        outfile.write("##gff-version 3\n")
        for input_gff in gffs:
            with gzip.open(input_gff, "rt") as ingff:
                for line in ingff:
                    line = line.rstrip()
                    if line.startswith("##FASTA"):
                        break
                    elif not line.startswith("#"):
                        writer.writerow(line.split("\t"))


def main():
    parser = argparse.ArgumentParser(
        description="Concatenate GFF files into a compressed TSV file."
    )
    parser.add_argument(
        "--gffs", nargs="+", help="List of input GFF files to concatenate"
    )
    parser.add_argument("--output", help="Output concatenated GFF file")

    args = parser.parse_args()

    concatenate_gff_files(args.gffs, args.output)


if __name__ == "__main__":
    main()
