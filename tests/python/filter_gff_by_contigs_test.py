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

import gzip
import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "bin"))

from filter_gff_by_contigs import contig_ids_from_fasta, filter_gff_by_contigs

DUMMY_SEQ = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"


def _cds(contig_id: str, feature_id: str, feature_type: str = "CDS") -> str:
    return (
        "\t".join(
            [
                contig_id,
                "pyrodigal_v3.6.3",
                feature_type,
                "1",
                "60",
                ".",
                "+",
                "0",
                f"ID={feature_id}",
            ]
        )
        + "\n"
    )


GFF_CONTENT = (
    "##gff-version 3\n"
    '# Sequence Data: seqnum=1;seqlen=120;seqhdr="ERZ12345_1"\n'
    + _cds("ERZ12345_1", "ERZ12345_1_1")  # keep
    + _cds("ERZ12345_1", "ERZ12345_1_2")  # keep
    + _cds("ERZ12345_2", "ERZ12345_2_1")  # exclude
    + _cds("ERZ12345_3", "ERZ12345_3_1_60_+")  # keep, FragGeneScanRS
    + _cds("ERZ12345_10", "ERZ12345_10_1")  # exclude, must not match ERZ12345_1
)


class TestContigIdsFromFasta(unittest.TestCase):
    """Tests for contig_ids_from_fasta."""

    def setUp(self):
        self.test_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_only_first_word_of_header_used(self):
        fasta = self.test_dir / "contigs.fasta"
        fasta.write_text(
            ">ERZ12345_1 len=120\n" + DUMMY_SEQ + ">ERZ12345_3\n" + DUMMY_SEQ
        )
        self.assertEqual(contig_ids_from_fasta(fasta), {"ERZ12345_1", "ERZ12345_3"})

    def test_compressed_fasta(self):
        fasta = self.test_dir / "contigs.fasta.gz"
        with gzip.open(fasta, "wt") as fh:
            fh.write(">ERZ12345_1\n" + DUMMY_SEQ)
        self.assertEqual(contig_ids_from_fasta(fasta), {"ERZ12345_1"})


class TestFilterGffByContigs(unittest.TestCase):
    """Tests for filter_gff_by_contigs."""

    def setUp(self):
        self.test_dir = Path(tempfile.mkdtemp())
        self.gff = self.test_dir / "assembly.gff.gz"
        with gzip.open(self.gff, "wt") as fh:
            fh.write(GFF_CONTENT)
        self.output = self.test_dir / "chunk.gff.gz"

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def _output_lines(self) -> list[str]:
        with gzip.open(self.output, "rt") as fh:
            return fh.read().splitlines()

    def test_keeps_only_the_contigs_in_the_chunk(self):
        cds_count = filter_gff_by_contigs(
            self.gff, {"ERZ12345_1", "ERZ12345_3"}, self.output
        )
        self.assertEqual(cds_count, 3)
        self.assertEqual(
            self._output_lines(),
            [
                "##gff-version 3",
                _cds("ERZ12345_1", "ERZ12345_1_1").rstrip("\n"),
                _cds("ERZ12345_1", "ERZ12345_1_2").rstrip("\n"),
                _cds("ERZ12345_3", "ERZ12345_3_1_60_+").rstrip("\n"),
            ],
        )

    def test_no_cds_for_the_contigs(self):
        cds_count = filter_gff_by_contigs(self.gff, {"ERZ12345_4"}, self.output)
        self.assertEqual(cds_count, 0)
        self.assertEqual(self._output_lines(), ["##gff-version 3"])

    def test_features_that_are_not_cds_are_not_counted(self):
        with gzip.open(self.gff, "wt") as fh:
            fh.write(
                "##gff-version 3\n"
                + _cds("ERZ12345_1", "ERZ12345_1_1", feature_type="gene")
            )
        cds_count = filter_gff_by_contigs(self.gff, {"ERZ12345_1"}, self.output)
        self.assertEqual(cds_count, 0)

    def test_uncompressed_input_and_output(self):
        gff = self.test_dir / "assembly.gff"
        gff.write_text(GFF_CONTENT)
        output = self.test_dir / "chunk.gff"
        cds_count = filter_gff_by_contigs(gff, {"ERZ12345_2"}, output)
        self.assertEqual(cds_count, 1)
        self.assertEqual(
            output.read_text(), "##gff-version 3\n" + _cds("ERZ12345_2", "ERZ12345_2_1")
        )


if __name__ == "__main__":
    unittest.main()
