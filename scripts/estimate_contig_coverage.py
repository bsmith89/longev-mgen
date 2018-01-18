#!/usr/bin/env python3

import sys
import pandas as pd

def calc_cvrg(df):
    total_count = df.depth.sum()
    assert len(df.length.unique()) == 1
    contig_len = df.length.iloc[0]
    hit_len = (df.depth > 0).sum()
    left = df.position.min()
    right = df.position.max()
    return pd.Series({'coverage': total_count / contig_len,
                      'hit_fraction': hit_len / contig_len,
                      'left': left,
                      'right': right,
                      'contiguous': (right - left) + 1 == hit_len})

def main():
    depth_path = sys.argv[1]
    contig_path = sys.argv[2]
    mgen_library_id = sys.argv[3]
    depth = pd.read_table(depth_path)
    contig = pd.read_table(contig_path, index_col='contig_id')

    data = depth.join(contig, on='contig_id').groupby('contig_id').apply(calc_cvrg)
    data['mgen_library_id'] = mgen_library_id
    data.left = data.left.astype(int)
    data.right = data.right.astype(int)
    data.contiguous = data.contiguous.astype(int)
    data.to_csv(sys.stdout, sep='\t')

if __name__ == '__main__':
    main()
