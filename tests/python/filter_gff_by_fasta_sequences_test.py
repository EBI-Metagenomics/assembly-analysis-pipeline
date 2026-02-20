#!/usr/bin/env python3

import os
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "bin"))

from filter_gff_by_fasta_sequences import (
    extract_fasta_sequence_ids,
    filter_gff_by_sequence_ids,
)


PRODIGAL_HEADER = (
    ">ERZ23524530_111_2 # 841 # 11490 # -1 # ID=12368_2;partial=01;start_type=Edge;\n"
)
FGSRS_HEADER_NEG = ">ERZ23524530_383157_129_194_-\n"
FGSRS_HEADER_POS = ">ERZ23524530_385119_2_249_+\n"
DUMMY_SEQ = "MTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQVNDVIKDLVKDLRDAGVHYNAK\n"

GFF_HEADER = "##gff-version 3\n# Test GFF file\n"
GFF_PRODIGAL = (
    "\t".join(
        [
            "ERZ23524530_111_2",
            "prodigal",
            "CDS",
            "100",
            "500",
            ".",
            "+",
            "0",
            "ID=ERZ23524530_111_2",
        ]
    )
    + "\n"
)
GFF_FGSRS_1 = (
    "\t".join(
        [
            "ERZ23524530_383157",
            "FragGeneScanRS",
            "CDS",
            "129",
            "194",
            ".",
            "-",
            ".",
            "ID=ERZ23524530_383157_129_194_-",
        ]
    )
    + "\n"
)
GFF_FGSRS_2 = (
    "\t".join(
        [
            "ERZ23524530_385119",
            "FragGeneScanRS",
            "CDS",
            "2",
            "249",
            ".",
            "+",
            ".",
            "ID=ERZ23524530_385119_2_249_+",
        ]
    )
    + "\n"
)
GFF_EXTRA = (
    "\t".join(
        [
            "ERZ23524530_999999",
            "prodigal",
            "CDS",
            "200",
            "600",
            ".",
            "+",
            "0",
            "ID=ERZ23524530_999999_1",
        ]
    )
    + "\n"
)


class TestExtractFastaSequenceIds(unittest.TestCase):
    """Tests for extract_fasta_sequence_ids."""

    def setUp(self):
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def _write_fasta(self, content: str) -> str:
        path = os.path.join(self.test_dir, "test.faa")
        with open(path, "w") as fh:
            fh.write(content)
        return path

    def test_prodigal_format(self):
        path = self._write_fasta(PRODIGAL_HEADER + DUMMY_SEQ)
        self.assertEqual(extract_fasta_sequence_ids(path), {"ERZ23524530_111_2"})

    def test_fraggenescanrs_negative_strand(self):
        path = self._write_fasta(FGSRS_HEADER_NEG + DUMMY_SEQ)
        self.assertEqual(
            extract_fasta_sequence_ids(path), {"ERZ23524530_383157_129_194_-"}
        )

    def test_fraggenescanrs_positive_strand(self):
        path = self._write_fasta(FGSRS_HEADER_POS + DUMMY_SEQ)
        self.assertEqual(
            extract_fasta_sequence_ids(path), {"ERZ23524530_385119_2_249_+"}
        )

    def test_mixed_formats(self):
        content = (
            PRODIGAL_HEADER
            + DUMMY_SEQ
            + FGSRS_HEADER_NEG
            + DUMMY_SEQ
            + FGSRS_HEADER_POS
            + DUMMY_SEQ
        )
        path = self._write_fasta(content)
        self.assertEqual(
            extract_fasta_sequence_ids(path),
            {
                "ERZ23524530_111_2",
                "ERZ23524530_383157_129_194_-",
                "ERZ23524530_385119_2_249_+",
            },
        )

    def test_invalid_header_raises(self):
        path = self._write_fasta(">UNRECOGNISED_FORMAT\n" + DUMMY_SEQ)
        with self.assertRaises(ValueError):
            extract_fasta_sequence_ids(path)

    def test_fraggenescanrs_too_few_parts_raises(self):
        path = self._write_fasta(">ERZ23524530_-\n" + DUMMY_SEQ)
        with self.assertRaises(ValueError):
            extract_fasta_sequence_ids(path)


class TestFilterGffBySequenceIds(unittest.TestCase):
    """Tests for filter_gff_by_sequence_ids."""

    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.gff_file = os.path.join(self.test_dir, "test.gff")
        self.output_file = os.path.join(self.test_dir, "filtered.gff")

        with open(self.gff_file, "w") as fh:
            fh.write(GFF_HEADER + GFF_PRODIGAL + GFF_FGSRS_1 + GFF_FGSRS_2 + GFF_EXTRA)

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def _read_output(self) -> str:
        with open(self.output_file) as fh:
            return fh.read()

    def test_returns_all_gff_gene_ids(self):
        gff_ids = filter_gff_by_sequence_ids(self.gff_file, set(), self.output_file)
        self.assertEqual(
            gff_ids,
            {
                "ERZ23524530_111_2",
                "ERZ23524530_383157_129_194_-",
                "ERZ23524530_385119_2_249_+",
                "ERZ23524530_999999_1",
            },
        )

    def test_filters_correct_sequences(self):
        keep = {
            "ERZ23524530_111_2",
            "ERZ23524530_383157_129_194_-",
            "ERZ23524530_385119_2_249_+",
        }
        filter_gff_by_sequence_ids(self.gff_file, keep, self.output_file)
        content = self._read_output()

        for seq_id in keep:
            self.assertIn(seq_id, content)
        self.assertNotIn("ERZ23524530_999999", content)

    def test_preserves_header_lines(self):
        filter_gff_by_sequence_ids(self.gff_file, set(), self.output_file)
        content = self._read_output()
        self.assertIn("##gff-version 3", content)
        self.assertIn("# Test GFF file", content)

    def test_malformed_gff_line_raises(self):
        malformed_gff = os.path.join(self.test_dir, "malformed.gff")
        with open(malformed_gff, "w") as fh:
            fh.write("only_three\tfields\there\n")
        with self.assertRaises(ValueError):
            filter_gff_by_sequence_ids(malformed_gff, set(), self.output_file)

    def test_missing_id_attribute_raises(self):
        no_id_gff = os.path.join(self.test_dir, "no_id.gff")
        with open(no_id_gff, "w") as fh:
            fh.write(
                "\t".join(
                    [
                        "ERZ23524530_111_2",
                        "prodigal",
                        "CDS",
                        "100",
                        "500",
                        ".",
                        "+",
                        "0",
                        "Name=foo",
                    ]
                )
                + "\n"
            )
        with self.assertRaises(ValueError):
            filter_gff_by_sequence_ids(no_id_gff, set(), self.output_file)


if __name__ == "__main__":
    unittest.main()
