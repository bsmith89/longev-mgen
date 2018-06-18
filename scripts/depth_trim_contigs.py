#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pandas as pd
from Bio.SeqIO import index, write
import numpy as np
from tqdm import tqdm

def depth_trim(data, thresh, window):
    begin = data.position.min()
    end = data.position.max()
    left_most = None
    right_most = None
    for i in range(begin, end):
        left = i - window / 2
        right = i + window / 2
        if data[(left < data.position) & (right > data.position)].depth.mean() > thresh:
            left_most = i
            break
    for i in range(end, begin, -1):
        left = i - window / 2
        right = i + window / 2
        if data[(left < data.position) & (right > data.position)].depth.sum() / window > thresh:
            right_most = i
            break
    if (left_most is None) or (right_most is None):
        return None
    else:
        return (left_most, right_most)


if __name__ == "__main__":
    seqs = index(sys.argv[1], 'fasta')
    data = pd.read_table(sys.argv[2], names=['contig_id', 'position', 'depth'],
                          index_col='contig_id')
    thresh = float(sys.argv[3])
    window_size = int(sys.argv[4])
    min_length = int(sys.argv[5])
    median_depth = data.depth.median()
    tally_seqs = 0
    tally_nucs = 0
    for contig_id in tqdm(seqs.keys()):
        trim = depth_trim(data.loc[contig_id], median_depth * thresh, window_size)
        if trim:
            left, right = trim
            if (right - left) > min_length:
                tally_seqs += 1
                tally_nucs += right - left
                write(seqs[contig_id][left:right], sys.stdout, 'fasta')
    print("Output {} sequences with {} positions.".format(tally_seqs, tally_nucs), file=sys.stderr)

