#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
USAGE: plot_nucmer_comparison.py <v1_v2.nucmer.coords.tsv> <v1.length.tsv> <v2.length.tsv> output.[png|pdf]

<v1_v2.nucmer.coords.tsv> : A processed version of Mummer4's show-coords -B format (see column_names)
<[v1|v2].length.tsv> : TSVs of contig lengths

"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import collections  as mc
from matplotlib.ticker import StrMethodFormatter
import sys

if __name__ == "__main__":
    column_names = [  'contig_id_2'         # [ 1]
                    , 'length_2'            # [ 3]
                    , 'contig_id_1'         # [ 6]
                    , 'start_2'             # [ 7] Start of alignment in the second variant's sequence
                    , 'end_2'               # [ 8] End of alignment in the second variant's sequence
                    , 'start_1'             # [ 9] Start of alignment in the first variant's sequence
                    , 'end_1'               # [10] End of alignment in the second variant's sequence
                    , 'aligned_identity'    # [11] Fraction of positions identical in alignment
                    , 'strand'              # [18] Is the second alignment inverted from the first?
                    , 'length_1'            # [19] 
                ]

    data = pd.read_table(sys.argv[1],
                        names=column_names)

    # Calculate various quantities
    data['alength_1'] = (data.end_1 - data.start_1).abs()
    data['alength_2'] = (data.end_2 - data.start_2).abs()
    data['coverage_1'] = data.alength_1 / data.length_1
    data['coverage_2'] = data.alength_2 / data.length_2


    # Make strand data something we can work with.
    data.strand = data.strand.map({'Plus': +1, 'Minus': -1}).fillna(+1)
    # If the inversion covers more than half of the contig, than it's not really an inversion, is it?


    # Change to python-style indexing
    data[['start_1', 'start_2']] -= 1

    # Add entries for unaligned contigs
    contigs_1 = (pd.read_table(sys.argv[2])
                    .rename(columns={'contig_id': 'contig_id_1', 'length': 'length_1'}))
    contigs_2 = (pd.read_table(sys.argv[3])
                    .rename(columns={'contig_id': 'contig_id_2', 'length': 'length_2'}))
    data = pd.merge(data, contigs_1, how='outer', on='contig_id_1', suffixes=('', '_1'))
    data = pd.merge(data, contigs_2, how='outer', on='contig_id_2', suffixes=('', '_2'))

    # Check data integrity and then transfer contig lengths into correct columns.
    assert data.dropna(subset=['length_1']).apply(lambda x: x.length_1 == x.length_1_1, axis=1).all()
    assert data.dropna(subset=['length_2']).apply(lambda x: x.length_2 == x.length_2_2, axis=1).all()
    data.length_1.fillna(data.length_1_1, inplace=True)
    data.length_2.fillna(data.length_2_2, inplace=True)
    data.drop(['length_1_1', 'length_2_2'], axis=1, inplace=True)

    # Everything except the contig_ids can be replaced with 0 when there's no alignment.
    _fillna_cols = ['start_1', 'end_1', 'start_2', 'end_2',
                    'aligned_identity', 'length_1', 'length_2'
                ]
    data[_fillna_cols] = data[_fillna_cols].fillna(0)

    # Sort alignments with long contigs first.
    data['combined_length'] = data.length_1 + data.length_2
    data.sort_values('combined_length', ascending=False, inplace=True)

    # Left hand side of alignment
    idx_right_1 = data[['contig_id_1', 'length_1']].drop_duplicates().set_index('contig_id_1').length_1.cumsum()
    idx_right_1.name = 'idx_right_1'
    idx_right_2 = data[['contig_id_2', 'length_2']].drop_duplicates().set_index('contig_id_2').length_2.cumsum()
    idx_right_2.name = 'idx_right_2'

    data = data.join(idx_right_1, on='contig_id_1')        # End of sequence_1 range
    data['idx_left_1'] = data.idx_right_1 - data.length_1  # Start of sequence_1 range
    data['idx_start_1'] = data.idx_left_1 + data.start_1   # Start of alignment_1 range
    data['idx_end_1'] = data.idx_left_1 + data.end_1       # End of alignment_1 range

    # Now the same for the query strand data
    data = data.join(idx_right_2, on='contig_id_2')        # End of sequence_2 range
    data['idx_left_2'] = data.idx_right_2 - data.length_2  # Start of sequence_2 range
    data['idx_start_2'] = data.idx_left_2 + data.start_2   # Start of alignment_2 range
    data['idx_end_2'] = data.idx_left_2 + data.end_2       # End of alignment_2 range

    d = data.copy()#[data.strand == -1]
    plt.style.use('dark_background')


    fig, ax = plt.subplots(figsize=(8, 8))
    line_table = list(zip(zip(d.idx_start_1, d.idx_start_2), zip(d.idx_end_1, d.idx_end_2)))
    lc = mc.LineCollection(line_table, color=d.strand.map({-1.: 'cyan', +1.: 'red'}).fillna('red'), lw=2)
    ax.add_collection(lc)
    #ax.axis('off')
    ax.tick_params(axis='x', rotation=-90)

    left_min = min(d.idx_left_1.min(), d.idx_left_2.min())
    right_max = max(d.idx_right_1.max(), d.idx_right_2.max())
    axis_length = right_max - left_min

    ax.set_xlim(left_min - 0.05 * axis_length,
                right_max + 0.05 * axis_length)
    ax.set_ylim(left_min - 0.05 * axis_length,
                right_max + 0.05 * axis_length)

    # Dummy artists for legend.
    art1 = plt.plot([], [], color='red', label='same')
    art2 = plt.plot([], [], color='cyan', label='inverted')
    plt.legend({'same': art1, 'inverted': art2}, loc='lower right')

    ax.xaxis.set_major_formatter(StrMethodFormatter('{x:0.1e}', ))
    ax.yaxis.set_major_formatter(StrMethodFormatter('{x:0.1e}', ))
    ax.set_aspect(aspect=1)
    ax.set_xlabel("Strain 1")
    ax.set_ylabel("Strain 2")

    fig.savefig(sys.argv[4])
