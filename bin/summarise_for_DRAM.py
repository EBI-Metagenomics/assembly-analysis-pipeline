import sys
import os
import pandas as pd
import glob
import gzip
import argparse

def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-s', "--ko_summaries", type=str, nargs='+', help="list of ko summaries")
    parser.add_argument('-k', "--ko_per_contigs", type=str, nargs='+', help="list of ko annotations")
    parser.add_argument('-i', "--interpro_summaries", type=str, nargs='+', help="list of interpro summaries")
    parser.add_argument('-d', "--dbcan_overviews", type=str, nargs='+', help="list of dbcan overview files")
    parser.add_argument('-p', "--prefix", type=str, help="file prefix")

    args = parser.parse_args(argv)

    return args

if __name__ == "__main__":
    args = parse_args(sys.argv[1:])

    analyses = glob.glob(os.path.join(args.input, "ERZ*"))
    assemblies_annotations = {}     # { assemblies_annotations[accessions] = set(contigs) }

    for assembly_analysis_path in analyses:

        # ----- Functional annotation files ----- #

        kegg_annotations = os.path.join(assembly_analysis_path, "assembly", "assembly_ko_per_contig.tsv.gz")
        kegg_description = os.path.join(assembly_analysis_path, "assembly", "assembly_ko_summary.tsv.gz")

        ips_annotations = os.path.join(assembly_analysis_path, "assembly", "assembly_interproscan.tsv")
        cazy_annotations = os.path.join(assembly_analysis_path, "assembly", "assembly_dbcan_overview.txt.gz")

        # TODO: the way I gather assembly accessions and study accession depend on the meta inherited from above
        assembly_accession = assembly_analysis_path.split('/')[-1]
        assemblies_annotations[assembly_accession] = []

        # ----- CAZy ----- #

        cazys = {}      # { cazys[contig] = cazy }

        if os.path.exists(cazy_annotations):
            with gzip.open(cazy_annotations, 'rt') as f:
                csv_reader = csv.reader(f, delimiter='\t')
                next(csv_reader)
                for line in csv_reader:
                    contig, _, hmmer, dbcan_sub, diamond, _ = line.strip().split('\t')
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
                    assemblies_annotations[assembly_accession].append(contig)

        # ----- kegg ----- #

        keggs = {}      # { keggs[contig] = kegg }
        keggs_description, contigs_keggs_description = {}, {}

        if os.path.exists(kegg_annotations) and os.path.exists(kegg_description):
            with gzip.open(kegg_annotations, 'rt') as f:
                csv_reader = csv.reader(f, delimiter='\t')
                for line in csv_reader:
                    ko, contig = line.strip().split('\t')
                    if contig not in keggs:
                        keggs[contig] = []
                    keggs[contig].append(ko)
                    assemblies_annotations[assembly_accession].append(contig)

            with gzip.open(kegg_description, 'rt') as f:
                csv_reader = csv.reader(f, delimiter='\t')
                for line in csv_reader:
                    _, ko, ko_description = line.strip().split('\t')
                    keggs_description[ko] = ko_description

            for contig in keggs:
                description_list = []
                if isinstance(keggs[contig], list):
                    for annotation in keggs[contig]:
                        description_list.append(keggs_description[annotation])
                    contigs_keggs_description[contig] = "; ".join(description_list)
                    keggs[contig] = "; ".join(keggs[contig])

        # ----- IPS ----- #

        pfams = {}      # { pfams[contig] = pfam }

        if os.path.exists(ips_annotations):
            with open(ips_annotations, 'r') as f: # convert to gzip once the file is compressed
                csv_reader = csv.reader(f, delimiter='\t')
                for line in csv_reader:
                    if "pfam" in line:
                        fields = line.strip().split('\t')
                        contig = fields[0]
                        ann_id = fields[4]
                        ann_desc = fields[5]
                        if contig not in pfams:
                            pfams[contig] = []
                        pfams[contig].append(ann_desc + " [" + ann_id + ']')

            for contig in pfams:
                pfams[contig] = "; ".join(pfams[contig])
                assemblies_annotations[assembly_accession].append(contig)

        # ----- Assemble tsv table ----- #

        table_header = ["", "fasta", "scaffold", "gene_position", "kegg_id", "kegg_hit", "pfam_hits", "cazy_id"]

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
                    keggs.get(contig, ""),                    # kegg_id
                    contigs_keggs_description.get(contig, "") # kegg_hit
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
