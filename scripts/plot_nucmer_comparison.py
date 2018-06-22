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
import matplotlib.patches as mp
from matplotlib.ticker import StrMethodFormatter
from matplotlib.backends.backend_pdf import PdfPages
from itertools import repeat
import sys

if __name__ == "__main__":
# FROM http://mummer.sourceforge.net/manual/#coords
# When run with the -B option, output format will consist of 21 tab-delimited
# columns. These are as follows: [1] query sequence ID [2] date of alignment
# [3] length of query sequence [4] alignment type [5] reference file [6]
# reference sequence ID [7] start of alignment in the query [8] end of
# alignment in the query [9] start of alignment in the reference [10] end of
# alignment in the reference [11] percent identity [12] percent similarity [13]
# length of alignment in the query [14] 0 for compatibility [15] 0 for
# compatibility [16] NULL for compatibility [17] 0 for compatibility [18]
# strand of the query [19] length of the reference sequence [20] 0 for
# compatibility [21] and 0 for compatibility.
    _column_names = [ 'contig_id_2'         # [ 1]
                    , '_1'
                    , 'length_2'            # [ 3]
                    , '_2'
                    , '_3'
                    , 'contig_id_1'         # [ 6]
                    , 'start_2'             # [ 7] Start of alignment in the second variant's sequence
                    , 'end_2'               # [ 8] End of alignment in the second variant's sequence
                    , 'start_1'             # [ 9] Start of alignment in the first variant's sequence
                    , 'end_1'               # [10] End of alignment in the second variant's sequence
                    , 'aligned_identity'    # [11] Fraction of positions identical in alignment
                    , '_4'
                    , '_5'
                    , '_6'
                    , '_7'
                    , '_8'
                    , '_9'
                    , 'strand'              # [18] Is the second alignment inverted from the first?
                    , 'length_1'            # [19] 
                    , '_10'
                    , '_11'
                ]
    column_names = [n for n in _column_names if not n.startswith('_')]

    data = pd.read_table(sys.argv[1], names=_column_names, usecols=column_names)
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
    data['sort_key'] = -data[['length_1', 'length_2']].max(1)
    data.sort_values('sort_key', inplace=True)

    # Left hand side of alignment
    idx_right_1 = data[['contig_id_1', 'length_1']].drop_duplicates().set_index('contig_id_1').length_1.cumsum()
    idx_right_1.name = 'idx_right_1'
    idx_right_2 = data[['contig_id_2', 'length_2']].drop_duplicates().set_index('contig_id_2').length_2.cumsum()
    idx_right_2.name = 'idx_right_2'

    data = data.join(idx_right_1, on='contig_id_1')        # End of sequence_1 range
    data['idx_left_1'] = data.idx_right_1 - data.length_1  # Start of sequence_1 range
    data['idx_start_1'] = data.idx_left_1 + data.start_1   # Start of alignment_1 range
    data['idx_end_1'] = data.idx_left_1 + data.end_1       # End of alignment_1 range
    data['idx_middle_1'] = (data.idx_start_1 + data.idx_end_1) / 2

    # Now the same for the query strand data
    data = data.join(idx_right_2, on='contig_id_2')        # End of sequence_2 range
    data['idx_left_2'] = data.idx_right_2 - data.length_2  # Start of sequence_2 range
    data['idx_start_2'] = data.idx_left_2 + data.start_2   # Start of alignment_2 range
    data['idx_end_2'] = data.idx_left_2 + data.end_2       # End of alignment_2 range
    data['idx_middle_2'] = (data.idx_start_2 + data.idx_end_2) / 2

    d = data.copy()#[data.strand == -1]


    # Constants
    color_inv = 'red'
    color_fwd = 'blue'
    d['color'] = d.strand.map({-1.: color_inv, +1.: color_fwd}).fillna(color_inv)
    d.dropna(subset=['contig_id_1', 'contig_id_2',
                     ], inplace=True)
    padding = 0.02

    assert sys.argv[4].rsplit('.', 1)[-1] == 'pdf', 'The output figure must be a PDF.'
    with PdfPages(sys.argv[4]) as pdf:
        # View 1 - "Dots"
        fig, ax = plt.subplots(figsize=(15, 15))
        line_table = list(zip(zip(d.idx_start_1, d.idx_start_2), zip(d.idx_end_1, d.idx_end_2)))
        lc = mc.LineCollection(line_table, color=d.color, lw=2)
        ax.add_collection(lc)
        ax.tick_params(axis='x', rotation=-90)

        left_min = min(d.idx_left_1.min(), d.idx_left_2.min())
        right_max = max(d.idx_right_1.max(), d.idx_right_2.max())
        axis_length = right_max - left_min

        ax.set_xlim(left_min - padding * axis_length,
                    right_max + padding * axis_length)
        ax.set_ylim(left_min - padding * axis_length,
                    right_max + padding * axis_length)

        # Dummy artists for legend.
        art_inv = plt.plot([], [], color=color_inv)
        art_fwd = plt.plot([], [], color=color_fwd)
        plt.legend({'same': art_fwd, 'inverted': art_inv}, loc='lower right')

        ax.xaxis.set_major_formatter(StrMethodFormatter('{x:0.1e}', ))
        ax.yaxis.set_major_formatter(StrMethodFormatter('{x:0.1e}', ))
        ax.set_aspect(aspect=1)
        ax.set_xlabel("Strain 1")
        ax.set_ylabel("Strain 2")
        pdf.savefig(fig)
        plt.close(fig)


        # View 2 - Lines
        fig, ax = plt.subplots(figsize=(15, 3))
        width_coef = 1.5e4  # Conversion from summed alength's to linewidths

        pcoords = (list(zip(zip(d.idx_start_1, repeat(0)), zip(d.idx_end_1, repeat(0)),
                        zip(d.idx_end_2, repeat(1)), zip(d.idx_start_2, repeat(1)))))
        pc = mc.PatchCollection([mp.Polygon(c, closed=True)
                                for c in pcoords],
                                color=d.strand.map({-1.: 'red', +1.: 'blue'}).fillna('red'),
                                alpha=0.4, lw=0)
        ax.add_collection(pc)

        left_min = min(d.idx_left_1.min(), d.idx_left_2.min())
        right_max = max(d.idx_right_1.max(), d.idx_right_2.max())
        axis_length = right_max - left_min

        ax.set_xlim(left_min - padding * axis_length,
                    right_max + padding * axis_length)
        ax.set_ylim(-0.1, 1.1)

        # Dummy artists for legend.
        art_inv = plt.plot([], [], color=color_inv)
        art_fwd = plt.plot([], [], color=color_fwd)
        plt.legend({'same': art_fwd, 'inverted': art_inv}, bbox_to_anchor=(1.12, 1))

        ax.xaxis.set_major_formatter(StrMethodFormatter('{x:0.1e}', ))
        ax.set_yticks([0, 1])
        ax.set_yticklabels(['Strain 1', 'Strain 2'])
        pdf.savefig(fig)
        plt.close(fig)
