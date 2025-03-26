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
import os
import argparse
import pandas as pd
import gzip
import csv

def get_accession(file_name):
    return file_name.split('_')[0]

def extract_kegg_annotations(kegg_anns):
    keggs = {}

    global assemblies_annotations

    for kegg_annotation in kegg_anns:
        analysis_accession = get_accession(kegg_annotation)
        assemblies_annotations[analysis_accession] = []
        with gzip.open(kegg_annotation, 'rt') as f:
            csv_reader = csv.reader(f, delimiter='\t')
            next(csv_reader)
            for ko, contig in csv_reader:
                if contig not in keggs:
                    keggs[contig] = []
                keggs[contig].append(ko)
                assemblies_annotations[analysis_accession].append(contig)

    return keggs

def extract_cazy_families(cazy_files):
    cazys = {}

    global assemblies_annotations

    for cazy_annotation in cazy_files:
        analysis_accession = get_accession(cazy_annotation)
        assemblies_annotations[analysis_accession] = []
        with gzip.open(cazy_annotation, 'rt') as f:
            csv_reader = csv.reader(f, delimiter='\t')
            next(csv_reader)
            for contig, _, hmmer, dbcan_sub, diamond, _ in csv_reader:
                hmmer_list = [s.split('(')[0] for s in hmmer.split('+')]         # format: GT9(116-281)+GT17(914-1160) ==> ['GT9', 'GT17']
                dbcan_sub_list = [s.split('_')[0] for s in dbcan_sub.split('+')] # format: GH23_e756+CBM50_e338 ==> ['GH23', 'CBM50']
                diamond_list = [s.split('_')[0] for s in diamond.split('+')]     # format: CBM50+GH23 ==> ['CBM50', 'GH23']
                consensus = set()
                for hmmer_acc in hmmer_list:
                    if hmmer_acc in dbcan_sub_list or hmmer_acc in diamond_list:
                        consensus.add(hmmer_acc)
                for dbcan in dbcan_sub_list:
                    if dbcan in diamond_list:
                        consensus.add(dbcan)
                cazys[contig] = "; ".join(consensus)
                assemblies_annotations[analysis_accession].append(contig)

    return cazys

def extract_pfam_annotations(ips_files):
    pfams = {}

    global assemblies_annotations

    for ips_annotation in ips_files:
        analysis_accession = get_accession(ips_annotation)
        assemblies_annotations[analysis_accession] = []
        with open(ips_annotation, 'r') as f: # convert to gzip once the file is compressed
            csv_reader = csv.reader(f, delimiter='\t')
            for line in csv_reader:
                if "Pfam" in line:
                    contig = line[0]
                    ann_id = line[4]
                    ann_desc = line[5]
                    if contig not in pfams:
                        pfams[contig] = []
                    pfams[contig].append(ann_desc + " [" + ann_id + ']')

        for contig in pfams:
            if isinstance(pfams[contig], list):
                pfams[contig] = "; ".join(pfams[contig])
                assemblies_annotations[analysis_accession].append(contig)

    return pfams

def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', "--interpro_summaries", type=str, nargs='+', help="list of interpro summaries")
    parser.add_argument('-k', "--ko_per_contigs", type=str, nargs='+', help="list of ko annotations")
    parser.add_argument('-d', "--dbcan_overviews", type=str, nargs='+', help="list of dbcan overview files")
    parser.add_argument('-p', "--prefix", type=str, help="file prefix")

    args = parser.parse_args(argv)

    return args

if __name__ == "__main__":
    args = parse_args(sys.argv[1:])

    assemblies_annotations = {}     # { assemblies_annotations[accessions] = set(contigs) }

    kegg_annotations = args.ko_per_contigs
    ips_annotations = args.interpro_summaries
    cazy_annotations = args.dbcan_overviews

    # ----- CAZy ----- #

    cazys = extract_cazy_families(cazy_annotations) # { cazys[contig] = cazy }
    
    # ----- kegg ----- #

    keggs = extract_kegg_annotations(kegg_annotations) # { keggs[contig] = kegg }

    # ----- IPS ----- #

    pfams = extract_pfam_annotations(ips_annotations)   # { pfams[contig] = pfam }

    # ----- Assemble tsv table ----- #

    table_header = ["", "fasta", "scaffold", "gene_position", "kegg_id", "pfam_hits", "cazy_id"]

    for annotation in assemblies_annotations:
        functional_summary = []
        all_contigs = assemblies_annotations[annotation]
        for contig in all_contigs:
            contig_annotations = []
            contig_annotations.extend([
                contig,
                annotation,                 # fasta
                annotation + "_" + contig,  # scaffold
                contig                      # gene_position
            ])
            contig_annotations.extend([
                keggs.get(contig, "") # kegg_id
            ])
            contig_annotations.append(
                pfams.get(contig, "") # pfam_hits
            )
            contig_annotations.append(
                cazys.get(contig, "") # cazy_id
            )

            functional_summary.append(contig_annotations)

    partial_matrix = pd.DataFrame(functional_summary, index=all_contigs, columns=table_header)

    try:
        output_matrix = pd.concat([output_matrix, partial_matrix], ignore_index=True)
    except NameError:
        output_matrix = partial_matrix

    output_matrix.to_csv(f"{args.prefix}_summary_for_DRAM.tsv", sep='\t', header=True, index=False)
