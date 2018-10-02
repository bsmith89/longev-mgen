#!/usr/bin/env python3

from Bio.SeqIO import index
from copy import deepcopy
import sys

if __name__ == "__main__":
    seq_index = index(sys.argv[1], 'fasta')
    for frag in sys.argv[2:]:
        seq_id, *indices = frag.rsplit(':', 1)
        rec = deepcopy(seq_index[seq_id])
        if indices:
            left, right = indices[0].split('-')
            left = int(left)
            right = int(right)
        else:
            left, right = 0, len(rec)
        if left > right:
            rec.seq = rec.seq[right:left].reverse_complement()
        else:
            rec.seq = rec.seq[left:right]
        print(f'>{rec.id}\n{rec.seq}')
