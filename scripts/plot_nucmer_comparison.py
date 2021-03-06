#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
USAGE: plot_nucmer_comparison.py <v1_v2.nucmer.coords.tsv> <v1.nlength.tsv> <v2.nlength.tsv> <alignment_length_thresh> output.[png|pdf]

<v1_v2.nucmer.coords.tsv> : A processed version of Mummer4's show-coords -B format (see column_names)
<[v1|v2].length.tsv> : TSVs of contig lengths

"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import collections as mc
import matplotlib.patches as mp
from matplotlib.ticker import StrMethodFormatter
from matplotlib.backends.backend_pdf import PdfPages
from itertools import repeat
import sys
import numpy as np

def flip_inverted_contigs(df, check=True):
    df = df.copy()
    assert df.contig_id_2.unique().shape[0] == 1
    if check:  # that the df should be flipped
        if not (df.alength_2 * df.strand).sum() < 0:
            return df
    df[['start_2', 'end_2']] = df[['end_2', 'start_2']]
    df.strand = df.strand * -1
    return df

def union_of_intervals(interval_list):
    interval_list = [(min(start, end), max(start, end)) for start, end in interval_list]
    out = []
    for start, end in sorted(interval_list):
        if out and out[-1][1] >= start - 1:
            out[-1][1] = max(out[-1][1], end)
        else:
            out.append([start, end])
    return out

def sum_aligned_length(interval_list):
    tally = 0
    for start, end in interval_list:
        tally += end - start
    return tally

def calc_unique_aligned_length(df):
    """Calculate total length aligned for each contig included in the data.

    *df* : data[['contig_id', 'start', 'end']].dropna()
    """
    out = (df.groupby('contig_id')
             .apply(lambda x: union_of_intervals(zip(x.start, x.end)))
             .apply(sum_aligned_length))
    return out

if __name__ == "__main__":
    data_path = sys.argv[1]
    length1_path = sys.argv[2]
    length2_path = sys.argv[3]
    cvrg1_path = sys.argv[4]
    cvrg2_path = sys.argv[5]
    lib1_path = sys.argv[6]
    lib2_path = sys.argv[7]
    athresh = int(sys.argv[8])
    out_path = sys.argv[9]

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

    data = pd.read_table(data_path, names=_column_names, usecols=column_names)
    # Calculate various quantities
    data['alength_1'] = (data.end_1 - data.start_1).abs()
    data['alength_2'] = (data.end_2 - data.start_2).abs()
    data['coverage_1'] = data.alength_1 / data.length_1
    data['coverage_2'] = data.alength_2 / data.length_2

    # Make strand data something we can work with.
    data.strand = data.strand.map({'Plus': +1, 'Minus': -1})
    # If the inversion covers more than half of the alignments,
    # than it's not really an inversion, is it?
    # Flip these alignments as though I took the reverse
    # complement of the entire contig.
    data = data.groupby(['contig_id_2']).apply(flip_inverted_contigs)

    # Change to python-style indexing
    # TODO: Is this right?  Do I need to move the end_i index too?
    data[['start_1', 'start_2']] -= 1

    # Filter alignments for length
    data = data[lambda x: (x.alength_1 > athresh) & (x.alength_2 > athresh)]

   # Add entries for unaligned contigs
    length_1 = (pd.read_table(length1_path)
                    .rename(columns={'contig_id': 'contig_id_1', 'length': 'length_1'}))
    length_2 = (pd.read_table(length2_path)
                    .rename(columns={'contig_id': 'contig_id_2', 'length': 'length_2'}))
    data = pd.merge(data, length_1, how='outer', on='contig_id_1', suffixes=('', '_1'))
    data = pd.merge(data, length_2, how='outer', on='contig_id_2', suffixes=('', '_2'))

    # Check data integrity and then transfer contig lengths into correct columns.
    assert data.dropna(subset=['length_1']).apply(lambda x: x.length_1 == x.length_1_1, axis=1).all()
    assert data.dropna(subset=['length_2']).apply(lambda x: x.length_2 == x.length_2_2, axis=1).all()
    data.length_1.fillna(data.length_1_1, inplace=True)
    data.length_2.fillna(data.length_2_2, inplace=True)
    data.drop(['length_1_1', 'length_2_2'], axis=1, inplace=True)



    # Unaligned contigs get zero-length.
    _fillna_cols = ['length_1', 'length_2', 'alength_1', 'alength_2']
    data[_fillna_cols] = data[_fillna_cols].fillna(0)

    # Calculate total aligned length for every match
    alength_1 = data.groupby(['contig_id_1'])['alength_1'].sum()
    alength_1.name = 'total_alength_1'
    alength_2 = data.groupby(['contig_id_2'])['alength_2'].sum()
    alength_2.name = 'total_alength_2'
    data = data.join(alength_1, on='contig_id_1').join(alength_2, on='contig_id_2')

    # TODO: Optimize layout.
    # Approach: ("pick teams") where each round the next best contig is ordered
    # at the end of the sequence.

    # Sort alignments relative to Strain 1
    data['length_longer'] = data[['length_1', 'length_2']].max(1)
    data['total_alength_longer'] = data[['total_alength_1', 'total_alength_2']].max(1)
    data['alength_longer'] = data[['alength_1', 'alength_2']].max(1)
    data.sort_values(['alength_1', 'alength_2', 'start_1'], ascending=[False, False, True],
                     inplace=True)

    # Replace values with 0 after sorting, for layouts.
    _fillna_cols = ['start_1', 'end_1', 'start_2', 'end_2']
    data[_fillna_cols] = data[_fillna_cols].fillna(0)

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

    # Calculate how good the arrangement is.
    layout_loss = np.sqrt(((data.idx_left_1 - data.idx_left_2)**2 + (data.idx_right_1 - data.idx_right_2)**2).sum())
    print("Log layout loss:", np.log(layout_loss), file=sys.stderr)

    # Calculate how much alignment there is for each genome
    total_length_1 = data[['contig_id_1', 'length_1']].drop_duplicates().length_1.sum()
    total_alength_1 = (calc_unique_aligned_length(data[['contig_id_1', 'start_1', 'end_1']]
                                                    .rename(columns=lambda x: x[:-2]))
                                                .sum())
    total_length_2 = data[['contig_id_2', 'length_2']].drop_duplicates().length_2.sum()
    total_alength_2 = (calc_unique_aligned_length(data[['contig_id_2', 'start_2', 'end_2']]
                                                    .rename(columns=lambda x: x[:-2]))
                                                .sum())
    print('{} of {} nucleotides ({:0.1%}) aligned in Strain 1'
              .format(int(total_alength_1), int(total_length_1),
                      total_alength_1 / total_length_1),
          file=sys.stderr)
    print('{} of {} nucleotides ({:0.1%}) aligned in Strain 2'
              .format(int(total_alength_2), int(total_length_2),
                      total_alength_2 / total_length_2),
          file=sys.stderr)


    # Plotting
    # Constants
    color_inv = 'red'
    color_fwd = 'blue'
    contig_palette = ['blue', 'orange']
    data['color'] = data.strand.map({-1.: color_inv, +1.: color_fwd})
    padding = 0.02
    tick_length = 0.05

    # TODO: Do we want to remove this plot entirely?
    # # View 1 - "Dots"
    # d = data.copy()
    # d.dropna(subset=['contig_id_1', 'contig_id_2',
    #                 ], inplace=True)
    #
    # fig, ax = plt.subplots(figsize=(15, 15))
    # line_table = list(zip(zip(d.idx_start_1, d.idx_start_2), zip(d.idx_end_1, d.idx_end_2)))
    # lc = mc.LineCollection(line_table, color=d.color, lw=2)
    # ax.add_collection(lc)
    # ax.tick_params(axis='x', rotation=-90)
    #
    # left_min = min(d.idx_left_1.min(), d.idx_left_2.min())
    # right_max = max(d.idx_right_1.max(), d.idx_right_2.max())
    # axis_length = right_max - left_min
    #
    # ax.set_xlim(left_min - padding * axis_length,
    #             right_max + padding * axis_length)
    # ax.set_ylim(left_min - padding * axis_length,
    #             right_max + padding * axis_length)
    #
    # # Dummy artists for legend.
    # art_inv, *_ = plt.plot([], [], color=color_inv)
    # art_fwd, *_ = plt.plot([], [], color=color_fwd)
    # plt.legend([art_fwd, art_inv], ['same', 'inverted'], loc='lower right')
    #
    # ax.xaxis.set_major_formatter(StrMethodFormatter('{x:0.1e}', ))
    # ax.yaxis.set_major_formatter(StrMethodFormatter('{x:0.1e}', ))
    # ax.set_aspect(aspect=1)
    # ax.set_xlabel("Strain 1")
    # ax.set_ylabel("Strain 2")
    # pdf.savefig(fig)
    # plt.close(fig)

    # View 2 - Lines
    fig, ax = plt.subplots(figsize=(15, 3))
    width_coef = 1.5e4  # Conversion from summed alength's to linewidths

    # Plot alignments
    da = data.dropna(subset=['contig_id_1', 'contig_id_2',
                            ])
    pcoords = (list(zip(zip(da.idx_start_1, repeat(1)), zip(da.idx_end_1, repeat(1)),
                        zip(da.idx_end_2, repeat(2)), zip(da.idx_start_2, repeat(2)))))
    pc = mc.PatchCollection([mp.Polygon(c, closed=True)
                            for c in pcoords],
                            color=da.color,
                            alpha=0.4, lw=0)
    ax.add_collection(pc)

    # Load depth data
    with open(lib1_path) as lib1_handle, open(lib2_path) as lib2_handle:
        lib1 = [line.strip() for line in lib1_handle]
        lib2 = [line.strip() for line in lib2_handle]

    cvrg1 = pd.read_table(cvrg1_path, index_col='contig_id')
    cvrg_1_1 = cvrg1[list(set(lib1) & set(cvrg1.columns))].sum(1)
    cvrg_1_1 = cvrg_1_1 / cvrg_1_1.median()
    cvrg_1_1.name = 'cvrg_1_1'

    cvrg_1_2 = cvrg1[list(set(lib2) & set(cvrg1.columns))].sum(1)
    cvrg_1_2 = cvrg_1_2 / cvrg_1_2.median()
    cvrg_1_2.name = 'cvrg_1_2'

    cvrg2 = pd.read_table(cvrg2_path, index_col='contig_id')
    cvrg_2_1 = cvrg2[list(set(lib1) & set(cvrg2.columns))].sum(1)
    cvrg_2_1 = cvrg_2_1 / cvrg_2_1.median()
    cvrg_2_1.name = 'cvrg_2_1'

    cvrg_2_2 = cvrg2[list(set(lib2) & set(cvrg2.columns))].sum(1)
    cvrg_2_2 = cvrg_2_2 / cvrg_2_2.median()
    cvrg_2_2.name = 'cvrg_2_2'

    data = (data.join(cvrg_1_1, on='contig_id_1')
                .join(cvrg_1_2, on='contig_id_1')
                .join(cvrg_2_1, on='contig_id_2')
                .join(cvrg_2_2, on='contig_id_2'))


    # Plot contigs
    d1 = data[['contig_id_1', 'idx_left_1', 'idx_right_1',
            'length_1', 'cvrg_1_1', 'cvrg_1_2']].dropna().drop_duplicates()
    d2 = data[['contig_id_2', 'idx_left_2', 'idx_right_2',
            'length_2', 'cvrg_2_1', 'cvrg_2_2']].dropna().drop_duplicates()

    poly_table_1_1 = list(zip(zip(d1.idx_left_1, 1 - np.sqrt(d1.cvrg_1_1) / 20),
                        zip(d1.idx_left_1, repeat(1)),
                        zip(d1.idx_right_1, repeat(1)),
                        zip(d1.idx_right_1, 1 - np.sqrt(d1.cvrg_1_1) / 20)
                        )
                    )
    pc_1_1 = mc.PolyCollection(poly_table_1_1, lw=0,
                            alpha=0.5, color='grey')
    ax.add_collection(pc_1_1)
    poly_table_1_2 = list(zip(zip(d1.idx_left_1, 1 - np.sqrt(d1.cvrg_1_2) / 20),
                        zip(d1.idx_left_1, repeat(1)),
                        zip(d1.idx_right_1, repeat(1)),
                        zip(d1.idx_right_1, 1 - np.sqrt(d1.cvrg_1_2) / 20)
                        )
                    )
    pc_1_2 = mc.PolyCollection(poly_table_1_2, lw=0,
                            alpha=1, color=contig_palette)
    ax.add_collection(pc_1_2)

    poly_table_2_2 = list(zip(zip(d2.idx_left_2, 2 + np.sqrt(d2.cvrg_2_2) / 20),
                        zip(d2.idx_left_2, repeat(2)),
                        zip(d2.idx_right_2, repeat(2)),
                        zip(d2.idx_right_2, 2 + np.sqrt(d2.cvrg_2_2) / 20)
                        )
                    )
    pc_2_2 = mc.PolyCollection(poly_table_2_2, lw=0,
                            alpha=0.5, color='grey')
    ax.add_collection(pc_2_2)
    poly_table_2_1 = list(zip(zip(d2.idx_left_2, 2 + np.sqrt(d2.cvrg_2_1) / 20),
                        zip(d2.idx_left_2, repeat(2)),
                        zip(d2.idx_right_2, repeat(2)),
                        zip(d2.idx_right_2, 2 + np.sqrt(d2.cvrg_2_1) / 20)
                        )
                    )
    pc_2_1 = mc.PolyCollection(poly_table_2_1, lw=0,
                            alpha=1, color=contig_palette)
    ax.add_collection(pc_2_1)

    left_min = min(data.idx_left_1.min(), data.idx_left_2.min())
    right_max = max(data.idx_right_1.max(), data.idx_right_2.max())
    axis_length = right_max - left_min

    ax.set_xlim(left_min - padding * axis_length,
                right_max + padding * axis_length)
    ax.set_ylim(0.9, 2.1)

    ax.xaxis.set_major_formatter(StrMethodFormatter('{x:0.1e}'))
    ax.set_yticks([1, 2])
    ax.set_yticklabels(['Strain 1', 'Strain 2'])
    fig.savefig(out_path)
