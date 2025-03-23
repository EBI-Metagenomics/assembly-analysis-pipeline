import sys
import os
import pandas as pd
import glob
import gzip
import argparse

def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', "--input", type=str, help="input folder")

    args = parser.parse_args(argv)

    return args

if __name__ == "__main__":
    args = parse_args(sys.argv[1:])

    analyses = glob.glob(os.path.join(args.input, "ERZ*"))
    MGYAs = {}     # { MGYAs[accessions] = set(contigs) }

    for assembly_analysis_path in analyses:

        # ----- Functional annotation files ----- #

        kegg_annotations = os.path.join(assembly_analysis_path, "assembly", "assembly_ko_per_contig.tsv.gz")
        kegg_description = os.path.join(assembly_analysis_path, "assembly", "assembly_ko_summary.tsv.gz")

        ips_annotations = os.path.join(assembly_analysis_path, "assembly", "assembly_interproscan.tsv")
        cazy_annotations = os.path.join(assembly_analysis_path, "assembly", "assembly_dbcan_overview.txt.gz")

        current_MGYA = assembly_analysis_path.split('/')[-1]
        MGYAs[current_MGYA] = []

        # ----- CAZy ----- #

        CAZys = {}      # { CAZys[contig] = CAZy }

        if os.path.exists(cazy_annotations):
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
                    CAZys[MGYA_contig] = "; ".join(consensus)
                    MGYAs[current_MGYA].append(MGYA_contig)

        # ----- KEGG ----- #

        KEGGs = {}      # { KEGGs[contig] = KEGG }
        KEGGs_description, contigs_KEGGs_description = {}, {}

        if os.path.exists(kegg_annotations) and os.path.exists(kegg_description):
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

        if os.path.exists(ips_annotations):
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

        table_header = ["", "fasta", "scaffold", "gene_position", "kegg_id", "kegg_hit", "pfam_hits", "cazy_id"]

        for MGYA in MGYAs:
            functional_summary = []
            all_contigs = MGYAs[MGYA]
            for contig in all_contigs:
                contig_annotations = []
                contig_annotations.extend([
                    contig,
                    MGYA + "_" + contig,  # fasta
                    MGYA,                 # scaffold
                    contig                # gene_position
                ])
                try:
                    contig_annotations.extend([
                        KEGGs[contig],                    # kegg_id
                        contigs_KEGGs_description[contig] # kegg_hit
                    ])
                except KeyError:
                    contig_annotations.extend(["", ""])
                try:
                    contig_annotations.append(
                        Pfams[contig] # pfam_hits
                    )
                except KeyError:
                    contig_annotations.append("")
                try:
                    contig_annotations.append(
                        CAZys[contig] # cazy_id
                    )
                except KeyError:
                    contig_annotations.append("")

                functional_summary.append(contig_annotations)

        partial_matrix = pd.DataFrame(functional_summary, index=all_contigs, columns=table_header)

        try:
            output_matrix = pd.concat([output_matrix, partial_matrix], ignore_index=True)
        except NameError:
            output_matrix = partial_matrix

    output_matrix.to_csv("summary_for_DRAM_split.tsv", sep='\t', header=True, index=False)
