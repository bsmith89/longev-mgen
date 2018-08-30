#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# USAGE: format_signalp_data.py [SIGNALP] [FASTA] > result.tsv

import sys
import pandas as pd
from Bio.SeqIO import index


def relative_pos_closest(string, pos, char):
    """Relative position of the closest occurence of a character in a string.

    """
    left = string[:pos]
    right = string[pos:]
    try:
        left_dist = pos - left.rindex(char)
    except ValueError:
        left_dist = None
    try:
        right_dist = right.index(char)
    except ValueError:
        right_dist = None
    if left_dist is None:
        return right_dist
    elif right_dist is None:
        return -left_dist
    elif left_dist >= right_dist:
        return right_dist
    else:
        return -left_dist

if __name__ == "__main__":
    data = pd.read_table(sys.argv[1], sep='\s+', comment='#',
                         names=['feature_id', '_1', '_2', '_3', 'position',
                                '_5', '_6', '_7', 'score', '_9', '_10', '_11'],
                         index_col='feature_id',
                         usecols=['feature_id', 'position', 'score'])
    seqindex = index(sys.argv[2], 'fasta')

    for feature_id, (_position, score) in data.iterrows():
        position = int(_position)
        assert position == _position
        position -= 1
        seq = seqindex[feature_id].seq
        c_pos = relative_pos_closest(str(seq), position, 'C')
        if c_pos is None:
            # Print an empty string if there is no 'C'.
            c_pos = ''

        print(feature_id, position, round(score, 3),
              # Where is the closest C relative to the cleavage site?
              c_pos,
              sep='\t')
