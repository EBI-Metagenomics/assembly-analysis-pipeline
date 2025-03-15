#!/usr/bin/env python

import sys
import csv
from Bio import SearchIO

def parse_hmmscan(input_handle, output_handle):
    writer = csv.writer(output_handle, delimiter='\t')
    for qresult in SearchIO.parse(input_handle, 'hmmer3-tab'):
        contig_id = qresult.id
        for hit in qresult.hits:
            writer.writerow([hit.id, contig_id, hit.description])

if __name__ == '__main__':
    """
    This script parses the output of the hmmscan - tblout and writes a subset as a TSV

    It reads from the standard input and writes to standard output.

    The output format is: KO_ID<tab>CONTIG_ID<tab>"KO_DESCRIPTION"
    """
    parse_hmmscan(sys.stdin, sys.stdout)
