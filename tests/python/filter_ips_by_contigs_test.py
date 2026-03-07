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
import io
import os
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "bin"))

from filter_ips_by_contigs import (
    _contig_id_from_protein_id,
    _contig_ids_from_fasta,
    _iter_filtered_rows,
    filter_ips,
)

DUMMY_SEQ = "MTEIAAMMVKELRESTGAGMMDCKNALSETNGDFDKAVQVNDVIKDLVKDLRDAGVHYNAK\n"

FASTA_CONTENT = ">ERZ12345_1 some description\n" + DUMMY_SEQ + ">ERZ12345_2\n" + DUMMY_SEQ


def _ips_row(protein_id: str) -> str:
    return (
        "\t".join(
            [
                protein_id,
                "abc123",
                "100",
                "200",
                "1e-10",
                "Pfam",
                "PF00001",
                "Domain",
                "-",
                "-",
                "-",
                "IPR000001",
                "Some domain",
                "GO:0001",
                "1.1.1.1",
            ]
        )
        + "\n"
    )


IPS_CONTENT = (
    _ips_row("ERZ12345_1_3")          # Pyrodigal — contig ERZ12345_1 — keep
    + _ips_row("ERZ12345_1_7")        # Pyrodigal — contig ERZ12345_1 — keep
    + _ips_row("ERZ12345_2_1")        # Pyrodigal — contig ERZ12345_2 — keep
    + _ips_row("ERZ12345_1_2_70_-")   # FGS       — contig ERZ12345_1 — keep
    + _ips_row("ERZ12345_2_589_657_-")# FGS       — contig ERZ12345_2 — keep
    + _ips_row("ERZ12345_3_1")        # Pyrodigal — contig ERZ12345_3 — exclude
    + _ips_row("ERZ12345_99_2_70_+")  # FGS       — contig ERZ12345_99 — exclude
)


CONTIG_IDS = {"ERZ12345_1", "ERZ12345_2"}


class TestContigIdFromProteinId(unittest.TestCase):
    """Tests for _contig_id_from_protein_id."""

    def test_pyrodigal_format(self):
        self.assertEqual(_contig_id_from_protein_id("ERZ12345_1_3", CONTIG_IDS), "ERZ12345_1")

    def test_fgs_format_negative_strand(self):
        self.assertEqual(_contig_id_from_protein_id("ERZ12345_1_2_70_-", CONTIG_IDS), "ERZ12345_1")

    def test_fgs_format_positive_strand(self):
        self.assertEqual(_contig_id_from_protein_id("ERZ12345_2_589_657_+", CONTIG_IDS), "ERZ12345_2")

    def test_unknown_contig_returns_none(self):
        self.assertIsNone(_contig_id_from_protein_id("ERZ12345_99_1", CONTIG_IDS))

    def test_fgs_unknown_contig_returns_none(self):
        self.assertIsNone(_contig_id_from_protein_id("ERZ12345_99_2_70_+", CONTIG_IDS))


class TestContigIdsFromFasta(unittest.TestCase):
    """Tests for _contig_ids_from_fasta."""

    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.fasta = os.path.join(self.test_dir, "contigs.fasta")

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def _write(self, content: str) -> str:
        with open(self.fasta, "w") as fh:
            fh.write(content)
        return self.fasta

    def test_returns_contig_ids(self):
        self.assertEqual(
            _contig_ids_from_fasta(self._write(FASTA_CONTENT)),
            {"ERZ12345_1", "ERZ12345_2"},
        )

    def test_only_first_token_of_header_used(self):
        self.assertEqual(
            _contig_ids_from_fasta(self._write(">ERZ12345_1 extra fields\n" + DUMMY_SEQ)),
            {"ERZ12345_1"},
        )

    def test_empty_fasta_exits(self):
        with self.assertRaises(SystemExit):
            _contig_ids_from_fasta(self._write(""))


class TestIterFilteredRows(unittest.TestCase):
    """Tests for _iter_filtered_rows."""

    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.ips = os.path.join(self.test_dir, "ips.tsv")

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def _write(self, content: str) -> str:
        with open(self.ips, "w") as fh:
            fh.write(content)
        return self.ips

    def test_filters_to_matching_contigs(self):
        rows = list(_iter_filtered_rows(self._write(IPS_CONTENT), CONTIG_IDS))
        protein_ids = [r[0] for r in rows]
        self.assertEqual(
            set(protein_ids),
            {"ERZ12345_1_3", "ERZ12345_1_7", "ERZ12345_2_1", "ERZ12345_1_2_70_-", "ERZ12345_2_589_657_-"},
        )

    def test_empty_contig_set_yields_nothing(self):
        self.assertEqual(list(_iter_filtered_rows(self._write(IPS_CONTENT), set())), [])

    def test_fgs_protein_matched(self):
        row = _ips_row("ERZ12345_1_2_70_-")
        rows = list(_iter_filtered_rows(self._write(row), {"ERZ12345_1"}))
        self.assertEqual(rows[0][0], "ERZ12345_1_2_70_-")

    def test_comment_and_blank_lines_skipped(self):
        rows = list(_iter_filtered_rows(
            self._write("# comment\n\n" + _ips_row("ERZ12345_1_1")),
            {"ERZ12345_1"},
        ))
        self.assertEqual(len(rows), 1)

    def test_raw_line_preserved(self):
        row = _ips_row("ERZ12345_1_1")
        rows = list(_iter_filtered_rows(self._write(row), {"ERZ12345_1"}))
        self.assertEqual(rows[0][1], row)


class TestFilterIps(unittest.TestCase):
    """Integration tests for filter_ips."""

    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.fasta = os.path.join(self.test_dir, "contigs.fasta")
        self.ips = os.path.join(self.test_dir, "ips.tsv")
        self.out = os.path.join(self.test_dir, "out.tsv.gz")

        with open(self.fasta, "w") as fh:
            fh.write(FASTA_CONTENT)
        with open(self.ips, "w") as fh:
            fh.write(IPS_CONTENT)

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def _run(self) -> str:
        captured = io.StringIO()
        sys.stdout = captured
        try:
            filter_ips(self.fasta, self.ips, self.out)
        finally:
            sys.stdout = sys.__stdout__
        return captured.getvalue()

    def test_filtered_tsv_contains_matching_rows(self):
        self._run()
        with gzip.open(self.out, "rt") as fh:
            content = fh.read()
        self.assertIn("ERZ12345_1_3", content)
        self.assertIn("ERZ12345_1_2_70_-", content)
        self.assertIn("ERZ12345_2_589_657_-", content)
        self.assertNotIn("ERZ12345_3_1", content)
        self.assertNotIn("ERZ12345_99_2_70_+", content)

    def test_stdout_contains_matched_protein_ids(self):
        stdout = self._run()
        self.assertEqual(
            set(stdout.splitlines()),
            {"ERZ12345_1_3", "ERZ12345_1_7", "ERZ12345_2_1", "ERZ12345_1_2_70_-", "ERZ12345_2_589_657_-"},
        )


if __name__ == "__main__":
    unittest.main()
