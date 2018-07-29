#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pandas as pd
from Bio.SeqIO import index, write
import numpy as np
from tqdm import tqdm
from copy import deepcopy
import argparse

from tqdm import TqdmSynchronisationWarning
import warnings
warnings.simplefilter("ignore", TqdmSynchronisationWarning)

# TODO: This definitely needs a unit test.
def corr_trim(corr, thresh, flank, window, offset=0):
    assert window % 2 == 0, "Window size must be an even number."
    if len(corr) < flank:
        return []
    if np.max(corr) < thresh:
        return []
    left_most = None
    right_most = None
    corr = np.asarray(corr)
    # Trim flanks in from the left.
    for i in range(0, len(corr) - flank + 1):
        left = i
        right = i + flank
        if corr[left:right].mean() > thresh:
            left_most = i
            break
    if left_most is None:
        # We didn't find any part of the sequence with enough corr.
        return []

    # Trim flanks in from the right.
    for i in range(len(corr), flank - 1, -1):
        left = i - flank
        right = i
        if corr[left:right].mean() > thresh:
            right_most = i
            break
    if right_most is None:
        # We didn't find any part of the sequence with enough corr.
        return []

    if (right_most - left_most) < window:
        # The sequence is too short to try and split.
        return [(left_most + offset, right_most + offset)]

    # Find internal breakpoint and recursively call corr_trim.
    break_at = None
    for i in range(left_most + window // 2, right_most - window // 2):
        left = i - window // 2
        right = i + window // 2
        if corr[left:right].mean() < thresh:
            # Found the first breakpoint scanning from the left.
            break_at = i
            break
    if break_at is None:
        # Found no internal breakpoints.
        return [(left_most + offset, right_most + offset)]
    else:
        # Found internal breakpoints recursively.
        # Concatenate and return.
        left_frags_trimmed = corr_trim(corr[left_most:break_at], thresh, flank, window, offset=left_most + offset)
        right_frags_trimmed = corr_trim(corr[break_at:right_most], thresh, flank, window, offset=break_at + offset)
        return left_frags_trimmed + right_frags_trimmed


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--corr-thresh", "-d", type=float,
                   default=0.75,
                   help="Correlation score at which to trim/split contigs [%(default)s]")
    p.add_argument("--flank-size", "-f", type=int, dest='flank_size',
                   default=100, help="Flank size [%(default)s]")
    p.add_argument("--window-size", "-w", type=int, dest='window_size',
                   default=100, help="Window size [%(default)s]")
    p.add_argument("--min-length", "-l", type=int, dest='min_length',
                   default=1000, help="Minimum contig length to output [%(default)s]")
    p.add_argument('seq_path', type=str, metavar="FASTA",
                   help="Sequences to be trimmed.")
    p.add_argument('corr_handle', type=argparse.FileType('r'), metavar="CORR",
                   help="Correlation table from calculate_per_position_stats.py")
    args = p.parse_args()

    seqs = index(args.seq_path, 'fasta')
    data = pd.read_table(args.corr_handle,
                         names=['contig_id', 'position', 'total_depth', 'cosine_similarity'])
    data.position = data.position - 1  # Convert to zero-indexed.
    tally_seqs = 0
    tally_nucs = 0
    for contig_id in tqdm(list(seqs.keys())):
        seq = seqs[contig_id].seq
        if len(seq) < args.min_length:
            continue
        d = data[data.contig_id == contig_id]
        if d.empty:
            print("\rWARNING: {} not found in corr data.".format(contig_id),
                  file=sys.stderr)
            continue
        corr = np.zeros(d.position.max() + 1)
        corr[d.position] = d.cosine_similarity
        for frag in corr_trim(corr, args.corr_thresh, args.flank_size, args.window_size):
            left, right = frag
            if (right - left) > args.min_length:
                tally_seqs += 1
                tally_nucs += right - left
                name = contig_id + '_{}_{}'.format(left, right)
                print('>{}\n{}'.format(name, seq[left:right]), file=sys.stdout)
    print("Output {} sequences with {} positions.".format(tally_seqs, tally_nucs), file=sys.stderr)
