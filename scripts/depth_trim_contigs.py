#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pandas as pd
from Bio.SeqIO import index, write
import numpy as np
from tqdm import tqdm
from copy import deepcopy
import argparse

def depth_trim(depth, thresh, flank, window, offset=0):
    assert window % 2 == 0, "Window size must be an even number."
    if len(depth) < flank:
        return []
    left_most = None
    right_most = None
    depth = np.asarray(depth)
    # Trim in from the left.
    for i in range(0, len(depth) - flank + 1):
        left = i
        right = i + flank
        if depth[left:right].mean() > thresh:
            left_most = i
            break
    # Trim in from the right.
    for i in range(len(depth), flank - 1, -1):
        left = i - flank
        right = i
        if depth[left:right].mean() > thresh:
            right_most = i
            break

    if (left_most is None) or (right_most is None):
        # We didn't find any part of the sequence with enough depth.
        return []
    elif (right_most - left_most) < window:
        # The sequence is too short to try and split.
        return [(left_most, right_most)]
    else:
        # Find internal breakpoint and recursively call depth_trim.
        break_at = None
        for i in range(left_most + window // 2, right_most - window // 2):
            left = i - window // 2
            right = i + window // 2
            if depth[left:right].mean() < thresh:
                break_at = i
        if break_at is None:
            # Found no internal breakpoints.
            return [(left_most + offset, right_most + offset)]
        else:
            # Found internal breakpoints recursively; concatenate and return.
            left_frags_trimmed = depth_trim(depth[left_most:break_at], thresh, flank, window, offset=left_most + offset)
            right_frags_trimmed = depth_trim(depth[break_at:right_most], thresh, flank, window, offset=break_at + offset)
            return left_frags_trimmed + right_frags_trimmed


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--depth-thresh", "-d", type=float, dest='thresh',
                   default=0.01,
                   help="Fraction of median depth at which to trim/split contigs [%(default)s]")
    p.add_argument("--flank-size", "-f", type=int, dest='flank_size',
                   default=100, help="Flank size [%(default)s]")
    p.add_argument("--window-size", "-w", type=int, dest='window_size',
                   default=100, help="Window size [%(default)s]")
    p.add_argument("--min-length", "-l", type=int, dest='min_length',
                   default=1000, help="Minimum contig length to output [%(default)s]")
    p.add_argument('--depth-out', type=argparse.FileType('w'),
                   dest='depth_out_handle', metavar="DEPTH_OUT",
                   help="Output trimmed depth file")
    p.add_argument('seq_path', type=str, metavar="FASTA",
                   help="Sequences to be trimmed.")
    p.add_argument('depth_handle', type=argparse.FileType('r'),
                   metavar="DEPTH", help="Depth table from `samtools depth`")
    args = p.parse_args()

    seqs = index(args.seq_path, 'fasta')
    data = pd.read_table(args.depth_handle, names=['contig_id', 'position', 'depth'])
    data.position = data.position - 1  # Convert to zero-indexed.
    median_depth = data.depth.median()
    print("Median depth determined to be", median_depth, file=sys.stderr)
    tally_seqs = 0
    tally_nucs = 0
    for contig_id in tqdm(list(seqs.keys())):
        d = data[data.contig_id == contig_id]
        if d.empty:
            print("\rWARNING: {} not found in depth data.".format(contig_id),
                  file=sys.stderr)
            continue
        depth = np.zeros(d.position.max() + 1)
        depth[d.position] = d.depth
        for frag in depth_trim(depth, median_depth * args.thresh, args.flank_size, args.window_size):
            left, right = frag
            if (right - left) > args.min_length:
                tally_seqs += 1
                tally_nucs += right - left
                seq = seqs[contig_id].seq[left:right]
                name = contig_id + '_{}_{}'.format(left, right)
                print('>{}\n{}'.format(name, seq), file=sys.stdout)
                if args.depth_out_handle:
                    (d[(d.position >= left) & (d.position < right)]
                        .to_csv(args.depth_out_handle, sep='\t',
                                header=False, index=False))
    print("Output {} sequences with {} positions.".format(tally_seqs, tally_nucs), file=sys.stderr)

