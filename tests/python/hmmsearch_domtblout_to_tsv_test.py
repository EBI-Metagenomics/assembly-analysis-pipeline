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

from collections import Counter
import csv
import gzip
import io
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bin"))

from hmmsearch_domtblout_to_tsv import parse_hmmsearch


class TestParseHmmsearch(unittest.TestCase):
    """Ensure all domain hits are parsed across concatenated gzip chunks."""

    def test_concatenated_domain_tables_preserve_protein_ko_pairs(self) -> None:
        chunks = [
            [("K00001", "ERZ12345_1_1"), ("K00001", "ERZ12345_1_1")],
            [],
            [
                ("K00001", "ERZ12345_1_1"),
                ("K00001", "ERZ12345_1_2"),
                ("K00002", "ERZ12345_1_1"),
            ],
        ]
        expected = Counter(pair for chunk in chunks for pair in chunk)

        with tempfile.TemporaryDirectory() as tmpdir:
            members = []
            for index, pairs in enumerate(chunks):
                # domtblout uses column 1 for protein and column 4 for KO.
                table = "# target name accession tlen query name accession qlen\n"
                for ko, protein in pairs:
                    table += (
                        f"{protein} - 300 {ko} - 100 1e-20 160.0 0.0 "
                        "1 2 1e-20 1e-20 80.0 0.0 1 100 1 100 1 100 0.99 description\n"
                    )
                # A chunk's footer must not prevent parsing the next chunk's hits.
                table += "\n# Program: hmmsearch\n# [ok]\n"
                member = Path(tmpdir) / f"chunk_{index}.domtbl.gz"
                member.write_bytes(gzip.compress(table.encode()))
                members.append(member)

            for ordered in (members, list(reversed(members))):
                with self.subTest(order=[member.name for member in ordered]):
                    combined = Path(tmpdir) / "combined.domtbl.gz"
                    with combined.open("wb") as handle:
                        for member in ordered:
                            handle.write(member.read_bytes())
                    output = io.StringIO()
                    parse_hmmsearch(str(combined), output)
                    # The parser preserves repeated domains and pairs across chunks.
                    rows = csv.reader(io.StringIO(output.getvalue()), delimiter="\t")
                    self.assertEqual(Counter(tuple(row) for row in rows), expected)


if __name__ == "__main__":
    unittest.main()
