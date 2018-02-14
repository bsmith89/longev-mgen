#!/usr/bin/env python3

import sys
from Bio.SeqIO import parse

def fragment_seq(seq, length, overlap):
    if len(seq) < length:
        yield seq
    else:
        next_output = seq[:length]
        seq = seq[length - overlap:]
        while len(seq) >= length:
            yield next_output
            next_output = seq[:length]
            seq = seq[length - overlap:]
        yield next_output + seq


if __name__ == "__main__":
    length = int(sys.argv[1])
    overlap = int(sys.argv[2])
    for rec in parse(sys.argv[3], 'fasta'):
        for i, subseq in enumerate(fragment_seq(rec.seq, length, overlap)):
            print(f'>{rec.id}_{i}')
            print(subseq)
