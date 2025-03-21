import sys
import os
import pandas as pd
import glob
import gzip
import argparse

input_folder = "/hps/nobackup/rdf/metagenomics/service-team/users/germanab/dram_tests/ERZ15803840/assembly"

def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', "--input", type=str, help="input folder")

    args = parser.parse_args(argv)

    return args

if __name__ == "__main__":
    # args = parse_args(sys.argv[1:])

    # ----- Functional annotation files ----- #

    kegg_annotations = os.path.join(input_folder, "assembly_ko_per_contig.tsv.gz")
    kegg_description = os.path.join(input_folder, "assembly_ko_summary.tsv.gz")

    ips_annotations = os.path.join(input_folder, "assembly_interproscan.tsv")
    cazy_annotations = os.path.join(input_folder, "assembly_dbcan_overview.txt.gz")

    genome_summary_file = "/hps/nobackup/rdf/metagenomics/service-team/users/germanab/nf-modules/modules/ebi-metagenomics/dram/distill/tests/fixtures/dram_dbs/genome_summary_form.tsv"
    
    MGYAs = {}     # { MGYAs[accessions] = set(contigs) }
    current_MGYA = "ERZ15803840"
    MGYAs[current_MGYA] = []

    # ----- CAZy ----- #

    CAZys = {}      # { CAZys[contig] = CAZy }

    genome_summary = {} # CAZY family descriptions from CAZY db
    with open(genome_summary_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            genome_summary[fields[0]] = fields[1]

    with gzip.open(cazy_annotations, 'rt') as f:
        next(f)
        for line in f:
            MGYA_contig, _, HMMER, dbCAN_sub, DIAMOND, _ = line.strip().split('\t')
            HMMER_list = [s.split('(')[0] for s in HMMER.split('+')]
            dbCAN_sub_list = [s.split('_')[0] for s in dbCAN_sub.split('+')]
            DIAMOND_list = [s.split('_')[0] for s in DIAMOND.split('+')]
            consensus = set()
            for hmmer in HMMER_list:
                if hmmer in dbCAN_sub_list or hmmer in DIAMOND_list:
                    consensus.add(hmmer)
            for dbcan in dbCAN_sub_list:
                if dbcan in DIAMOND_list:
                    consensus.add(dbcan)
            CAZys[MGYA_contig] = ';'.join(consensus)
            MGYAs[current_MGYA].append(MGYA_contig)

    # for CDS in CAZy_summary:
    #     CAZYannotations = CAZy_summary[CDS].split(',')
    #     finalCAZYstring = ""
    #     for CAZy in CAZYannotations:
    #         try:
    #             notation = genome_summary[CAZy]
    #         except KeyError:
    #             notation = "None"
    #         finalCAZYstring += notation  + " [" + CAZy + ']'
    #         if len(CAZYannotations) > 1 and CAZYannotations.index(CAZy) != len(CAZYannotations)-1:
    #             finalCAZYstring += "; "
    #     CAZys[CDS] = finalCAZYstring

    # ----- KEGG ----- #
    
    KEGGs = {}      # { KEGGs[contig] = KEGG }
    KEGGs_description, contigs_KEGGs_description = {}, {}

    with gzip.open(kegg_annotations, 'rt') as f:
        for line in f:
            KO, MGYA_contig = line.strip().split('\t')
            if MGYA_contig not in KEGGs:
                KEGGs[MGYA_contig] = []
            KEGGs[MGYA_contig].append(KO)
            MGYAs[current_MGYA].append(MGYA_contig)

    with gzip.open(kegg_description, 'rt') as f:
        for line in f:
            _, KO, KO_description = line.strip().split('\t')
            KEGGs_description[KO] = KO_description

    for MGYA_contig in KEGGs:
        description_list = []
        if isinstance(KEGGs[MGYA_contig], list):
            for annotation in KEGGs[MGYA_contig]:
                description_list.append(KEGGs_description[annotation])
            contigs_KEGGs_description[MGYA_contig] = "; ".join(description_list)
            KEGGs[MGYA_contig] = "; ".join(KEGGs[MGYA_contig])

    # ----- IPS ----- #

    Pfams = {}      # { Pfams[contig] = Pfam }

    with open(ips_annotations, 'r') as f:
        for line in f:
            if "Pfam" in line:
                fields = line.strip().split('\t')
                MGYA_contig = fields[0]
                annID = fields[4]
                annDesc = fields[5]
                if MGYA_contig not in Pfams:
                    Pfams[MGYA_contig] = []
                Pfams[MGYA_contig].append(annDesc + " [" + annID + ']')
    
    for MGYA_contig in Pfams:
        Pfams[MGYA_contig] = "; ".join(Pfams[MGYA_contig])
        MGYAs[current_MGYA].append(MGYA_contig)
        
    # ----- Assemble tsv table ----- #

    table_header = ["", "fasta", "scaffold", "gene_position", "kegg_id", "kegg_hit", "pfam_hits", "cazy_hits"]
    functional_summary = []
    
    for MGYA in MGYAs:
        all_contigs = MGYAs[MGYA]
        for contig in all_contigs:
            rowList = []
            rowList.extend([
                contig,
                MGYA + "_" + contig,  # fasta
                MGYA,                 # scaffold
                contig                # gene_position
            ])
            try:
                rowList.extend([
                    KEGGs[contig],                    # kegg_id
                    contigs_KEGGs_description[contig] # kegg_hit
                ])
            except KeyError:
                rowList.extend(["", ""])
            try:
                rowList.append(
                    Pfams[contig] # pfam_hits
                )
            except KeyError:
                rowList.append("")
            try:
                rowList.append(
                    CAZys[contig] # kegg_hit
                )
            except KeyError:
                rowList.append("")

            functional_summary.append(rowList)

    # add mechanism to join more tbales together

    output_matrix = pd.DataFrame(functional_summary, index=all_contigs, columns=table_header)
    output_matrix.to_csv("summary_for_DRAM.tsv", sep='\t', header=True, index=False)
