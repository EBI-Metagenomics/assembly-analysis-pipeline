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
import re
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "bin"))

from summarise_for_dram import (
    extract_cazy_families_from_dbcan,
    extract_kegg_orthologs_annotations,
    extract_pfam_annotations,
)

# Regex that DRAM 1.3.5 uses to read the pfam_hits column
# https://github.com/WrightonLabCSU/DRAM/blob/91a902480c49adc52e991d898f73cf48ea3df174/mag_annotator/utils.py#L97-L99
DRAM_PFAM_REGEX = r"\[PF\d\d\d\d\d.\d*\]"

KO_PER_CONTIG_ROWS = [
    ["contig_id", "ko"],
    ["ERZ12345_1", "K00001", "K00002", "K00003"],
]

# fmt: off
DBCAN_OVERVIEW_ROWS = [
    ["Gene ID", "EC#", "dbCAN_hmm", "dbCAN_sub", "DIAMOND", "#ofTools", "Recommend Results"],
    # Three proteins of contig ERZ12345_1 have a consensus family
    ["ERZ12345_1_1", "-", "GT4(1-200)", "GT4_e10", "GT4", "3", "GT4"],
    ["ERZ12345_1_2", "-", "GT2(1-150)+CE4(160-300)", "GT2_e5+CE4_e1", "GT2", "3", "GT2+CE4"],
    ["ERZ12345_1_3", "-", "GH31(10-400)", "-", "GH31", "2", "GH31"],
    # Only one tool finds these families, so there is no consensus
    ["ERZ12345_1_4", "-", "GT51(1-90)", "-", "-", "1", "GT51"],
    ["ERZ12345_2_1", "-", "-", "-", "CBM50", "1", "CBM50"],
]
# fmt: on

# The script reads only the first 6 columns of the InterProScan TSV
IPS_ROWS = [
    ["ERZ12345_1_1", "abc123", "300", "Pfam", "PF00005", "ABC transporter"],
    ["ERZ12345_1_2", "abc123", "300", "Pfam", "PF00001", "7 transmembrane receptor"],
    ["ERZ12345_1_2", "abc123", "300", "PANTHER", "PTHR00001", "Not a Pfam entry"],
]


def _write_tsv_gz(test: unittest.TestCase, rows: list[list[str]]) -> str:
    """Write the rows to a gzip TSV file that is removed after the test."""
    tmp_dir = tempfile.TemporaryDirectory()
    test.addCleanup(tmp_dir.cleanup)
    path = os.path.join(tmp_dir.name, "input.tsv.gz")
    with gzip.open(path, "wt") as f:
        f.writelines("\t".join(row) + "\n" for row in rows)
    return path


class TestExtractKeggOrthologsAnnotations(unittest.TestCase):
    def test_kos_are_comma_separated_without_spaces(self):
        # DRAM module coverage splits on "," and does not strip the values
        keggs = extract_kegg_orthologs_annotations(
            _write_tsv_gz(self, KO_PER_CONTIG_ROWS)
        )
        self.assertEqual(keggs, {"ERZ12345_1": "K00001,K00002,K00003"})


class TestExtractCazyFamiliesFromDbcan(unittest.TestCase):
    def test_consensus_families_of_all_proteins_are_kept(self):
        # Before the fix, only the families of the last protein were kept
        cazys = extract_cazy_families_from_dbcan(
            _write_tsv_gz(self, DBCAN_OVERVIEW_ROWS)
        )
        self.assertEqual(cazys, {"ERZ12345_1": "CE4; GH31; GT2; GT4"})


class TestExtractPfamAnnotations(unittest.TestCase):
    def setUp(self):
        self.pfams = extract_pfam_annotations(_write_tsv_gz(self, IPS_ROWS))

    def test_only_pfam_with_the_placeholder_version(self):
        self.assertEqual(
            self.pfams,
            {
                "ERZ12345_1": "; ".join(
                    [
                        "ABC transporter [PF00005.1]",
                        "7 transmembrane receptor [PF00001.1]",
                    ]
                )
            },
        )

    def test_dram_regex_reads_all_accessions(self):
        # Without the version, the DRAM regex does not match the accession
        self.assertEqual(re.findall(DRAM_PFAM_REGEX, "ABC transporter [PF00005]"), [])
        dram_accessions = [
            match[1:-1].split(".")[0]
            for match in re.findall(DRAM_PFAM_REGEX, self.pfams["ERZ12345_1"])
        ]
        self.assertEqual(dram_accessions, ["PF00005", "PF00001"])


if __name__ == "__main__":
    unittest.main()
