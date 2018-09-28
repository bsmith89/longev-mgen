#!/usr/bin/env python3

from Bio.SeqIO import index, write
import sys

if __name__ == "__main__":
    rec_index = index(sys.argv[1], 'fasta')
    length_map = {k: len(rec_index[k].seq.ungap('-')) for k in rec_index}
    max_length = max(length_map.values())
    for k, length in length_map.items():
        if length >= max_length / 2:
            write(rec_index[k], sys.stdout, 'fasta')

