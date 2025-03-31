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

import argparse
import pandas as pd
import gzip
import csv
import Path
from collections import defaultdict


def extract_kegg_orthologs_annotations(kegg_summary_tsv: Path) -> dict:
    """Extract and gather the KO annotations per contig"""
    keggs = defaultdict(list())
    with gzip.open(kegg_summary_tsv, "rt") as f:
        csv_reader = csv.reader(f, delimiter="\t")
        next(csv_reader)
        for ko, contig in csv_reader:
            keggs[contig].append(ko)

    # There are no KOs per contig that are lists, right?
    for contig in keggs:
        if isinstance(keggs[contig], list):
            keggs[contig] = "; ".join(keggs[contig])

    return keggs


def extract_cazy_families_from_dbcan(dbcan_overview_tsv: Path) -> dict:
    """Extract and gather the CAZy families per contig from the run_DBCan overview file"""
    cazys = defaultdict(list())

    with gzip.open(dbcan_overview_tsv, "rt") as f:
        csv_reader = csv.reader(f, delimiter="\t")
        next(csv_reader)
        for contig, _, hmmer, dbcan_sub, diamond, _ in csv_reader:
            hmmer_list = [
                s.split("(")[0] for s in hmmer.split("+")
            ]  # format: GT9(116-281)+GT17(914-1160) ==> ['GT9', 'GT17']
            dbcan_sub_list = [
                s.split("_")[0] for s in dbcan_sub.split("+")
            ]  # format: GH23_e756+CBM50_e338 ==> ['GH23', 'CBM50']
            diamond_list = [
                s.split("_")[0] for s in diamond.split("+")
            ]  # format: CBM50+GH23 ==> ['CBM50', 'GH23']
            consensus = set()
            for hmmer_acc in hmmer_list:
                if hmmer_acc in dbcan_sub_list or hmmer_acc in diamond_list:
                    if not hmmer_acc == "-":
                        consensus.add(hmmer_acc)
            for dbcan in dbcan_sub_list:
                if dbcan in diamond_list:
                    if not dbcan == "-":
                        consensus.add(dbcan)

            if len(consensus) > 0:
                cazys[contig] = "; ".join(consensus)

    return cazys


def extract_pfam_annotations(interproscan_summary_tsv: Path) -> dict:
    pfams = defaultdict
    with gzip.open() as file_handler:
        csv_reader = csv.reader(file_handler, delimiter="\t")
        for line in csv_reader:
            contig, _, _, analysis, signature_accession, signature_description, *_ = (
                line
            )
            if analysis == "Pfam":
                pfams[contig].append(f"{signature_description} [{signature_accession}]")

    # There are no pfams per contig that are lists, right?
    for contig in pfams:
        if isinstance(pfams[contig], list):
            pfams[contig] = "; ".join(pfams[contig])

    return pfams


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-f",
        "--fasta",
        type=str,
        help="The fasta file with the contigs used for the analysis.",
    )
    parser.add_argument(
        "-i",
        "--interproscan",
        type=str,
        help="InterProScan TSV output (not the summary)",
    )
    parser.add_argument("-k", "--ko-per-contig", type=str, help="KO per contig TSV")
    parser.add_argument(
        "-d",
        "--dbcan-overview",
        type=str,
        help="run_DBCan overview output TSV file",
    )
    parser.add_argument("-p", "--prefix", type=str, help="TSV output prefix")

    args = parser.parse_args()

    assemblies_annotations = {}  # { assemblies_annotations[accessions] = set(contigs) }

    kegg_annotations = args.ko_per_contigs
    ips_annotations = args.interpro_summaries
    cazy_annotations = args.dbcan_overviews

    # ----- CAZy ----- #

    cazys = extract_cazy_families_from_dbcan(
        cazy_annotations
    )  # { cazys[contig] = cazy }

    # ----- KEGG ----- #

    keggs = extract_kegg_orthologs_annotations(
        kegg_annotations
    )  # { keggs[contig] = kegg }

    # ----- IPS ----- #

    pfams = extract_pfam_annotations(ips_annotations)  # { pfams[contig] = pfam }

    # ----- Assemble tsv table ----- #

    TABLE_HEADER = [
        "",
        "fasta",
        "scaffold",
        "gene_position",
        "kegg_id",
        "pfam_hits",
        "cazy_id",
    ]

    # Collect the contig names
    contigs = set()
    with open(args.fasta) as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                contigs.add(line.strip().replace(">", ""))

    functional_summary = []

    for contig in contigs:
        contig_annotations = []
        contig_annotations.extend(
            [
                contig,
                args.fasta,
                args.fasta + "_" + contig,  # scaffold
                contig,  # gene_position
            ]
        )
        contig_annotations.extend(
            [
                keggs.get(contig, "")  # kegg_id
            ]
        )
        contig_annotations.append(
            pfams.get(contig, "")  # pfam_hits
        )
        contig_annotations.append(
            cazys.get(contig, "")  # cazy_id
        )

        functional_summary.append(contig_annotations)

    partial_matrix = pd.DataFrame(
        functional_summary, index=contigs, columns=TABLE_HEADER
    )

    try:
        output_matrix = pd.concat([output_matrix, partial_matrix], ignore_index=True)
    except NameError:
        output_matrix = partial_matrix

    output_matrix.to_csv(
        f"{args.prefix}_summary_for_DRAM.tsv", sep="\t", header=True, index=False
    )
