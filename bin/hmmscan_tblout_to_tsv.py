#!/usr/bin/env python

import sys
from Bio import SearchIO

def parse_hmmscan(input_handle, output_handle):
    for qresult in SearchIO.parse(input_handle, 'hmmer3-tab'):
        for hit in qresult.hits:
            # This prints the hit id (KO id), the query id (the contig id) and the KO description
            output_handle.write(f"{hit.id}\t{hit.name}\t\"{hit.description}\"")

if __name__ == '__main__':
    """
    This script parses the output of the hmmscan - tblout and writes a subset as a TSV

    It reads from the standard input and writes to standard output.

    The output format is: KO_ID<tab>CONTIG_ID<tab>"KO_DESCRIPTION"
    """
    parse_hmmscan(sys.stdin, sys.stdout)
