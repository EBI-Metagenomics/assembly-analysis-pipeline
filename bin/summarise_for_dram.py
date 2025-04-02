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
from pathlib import Path
from collections import defaultdict
import functools


@functools.cache
def get_contig_id(fasta_header: str) -> str:
    """Proteins annotated with this pipeline have the following format:
    <assembly>_<contig_id>_<protein_id>
    """
    # TODO: add some validation here
    assembly, contig_id, *_ = fasta_header.split("_")
    return f"{assembly}_{contig_id}"


def extract_kegg_orthologs_annotations(aggregated_kos_per_contig: Path) -> dict:
    """Extract and gather the KO annotations per contig"""
    keggs = {}
    with gzip.open(aggregated_kos_per_contig, "rt") as f:
        csv_reader = csv.reader(f, delimiter="\t")
        next(csv_reader)
        # Example
        # contig_id       ko
        # ERZ1049444_314477_1     K00083  K17818  K22132  K13980
        for contig_id, ko_list in csv_reader:
            keggs[contig_id] = "; ".join( ko_list.strip().split(" ") )
    return keggs


def extract_cazy_families_from_dbcan(dbcan_overview_tsv: Path) -> dict:
    """Extract and gather the CAZy families per contig from the run_DBCan overview file"""
    cazys = defaultdict(list)

    with gzip.open(dbcan_overview_tsv, "rt") as f:
        csv_reader = csv.reader(f, delimiter="\t")
        next(csv_reader)
        for protein_id, _, hmmer, dbcan_sub, diamond, _ in csv_reader:
            contig = get_contig_id(protein_id)
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
    pfams = defaultdict(list)
    with gzip.open(interproscan_summary_tsv, "rt") as file_handler:
        csv_reader = csv.reader(file_handler, delimiter="\t")
        for line in csv_reader:
            # [
            #   'ERZ12345_4_1',
            #   _,
            #   _,
            #   'Pfam',
            #   'PF12574',
            #   '120 KDa Rickettsia surface antigen',
            #   ...
            (
                protein_id,
                _,
                _,
                analysis,
                signature_accession,
                signature_description,
                *_
            ) = line
            contig = get_contig_id(protein_id)
            if analysis.strip() == "Pfam":
                pfams[contig].append(f"{signature_description} [{signature_accession}]")

    for contig in pfams:
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
        help="The faa with the proteins of the assembly.",
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

    aggregated_kos_per_contig = args.ko_per_contig
    ips_annotations = args.interproscan
    dbcan_overview = args.dbcan_overview

    # ----- CAZy ----- #

    cazys = extract_cazy_families_from_dbcan(dbcan_overview)  # { cazys[contig] = cazy }

    # ----- KEGG ----- #

    keggs = extract_kegg_orthologs_annotations(
        aggregated_kos_per_contig
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

    # Collect the contig names from the proteins
    contigs = set()
    with gzip.open(args.fasta, "rt") as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                # The proteins have the <prefix>_<id>_[pyrodigal,fgsrs] at this point
                # So the contigs are named <prefix>_<incremental id>.. the prefix is
                # usually the assembly ERZ accession
                assembly, contig_id, *_ = line.strip().replace(">", "").split("_")
                contigs.add(f"{assembly}_{contig_id}")

    functional_summary = []

    for contig in contigs:
        contig_annotations = []
        contig_annotations.extend(
            [
                contig,
                args.prefix,  # assembly
                args.prefix + "_" + contig,  # scaffold
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

    partial_matrix = pd.DataFrame(functional_summary, columns=TABLE_HEADER)

    try:
        output_matrix = pd.concat([output_matrix, partial_matrix], ignore_index=True)
    except NameError:
        output_matrix = partial_matrix

    output_matrix.to_csv(
        f"{args.prefix}_summary_for_dram.tsv", sep="\t", header=True, index=False
    )
