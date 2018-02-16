#!/usr/bin/env python3
"""
{input.script} {input.bins} {input.merge_list} > {output}
"""

import pandas as pd
import sys

if __name__ == "__main__":
    bins_path = sys.argv[1]
    merge_path = sys.argv[2]
    min_complete_delta = float(sys.argv[3])
    max_contam_delta = float(sys.argv[4])

    bins = pd.read_table(bins_path, index_col='contig_id')
    merge = (pd.read_table(merge_path, skiprows=1,
                           names=['bin_A', 'bin_B',
                                  'complete_A', 'contam_A',
                                  'complete_B', 'contam_B',
                                  'complete_delta', 'contam_delta',
                                  'improvement_index',
                                  'complete_merge', 'contam_merge'])
                .drop(['complete_delta', 'contam_delta', 'improvement_index'],
                      axis='columns')
            )
    merge_flip = pd.DataFrame({'bin_A': merge.bin_B,
                               'bin_B': merge.bin_A,
                               'complete_A': merge.complete_B,
                               'contam_A': merge.contam_B,
                               'complete_B': merge.complete_A,
                               'contam_B': merge.contam_A,
                               'complete_merge': merge.complete_merge,
                               'contam_merge': merge.contam_merge})
    merge = pd.concat([merge, merge_flip])
    merge = merge[lambda x: x.complete_A > x.complete_B]
    merge['complete_delta'] = merge.complete_merge - merge.complete_A
    merge['contam_delta'] = merge.contam_merge - merge.contam_A
    merge['improvement_index'] = merge.complete_delta - merge.contam_delta
    merge = (merge.sort_values('improvement_index', ascending=False)
                  .groupby('bin_A').first()
                  .reset_index())
    merge = merge[(merge.complete_delta > min_complete_delta) &
                  (merge.contam_delta < max_contam_delta)]
    merge = merge[~merge.bin_B.isin(merge.bin_A)]
    assert not (set(merge.bin_A) & set(merge.bin_B))

    merge['metabin_name'] = merge.apply(lambda x: 'merge_{}_{}'.format(x.bin_A, x.bin_B), axis='columns')

    rename_map = {}
    for ix, row in merge.iterrows():
        rename_map[row.bin_A] = row.metabin_name
        rename_map[row.bin_B] = row.metabin_name

    bins.bin_id = bins.bin_id.replace(rename_map)

    bins.to_csv(sys.stdout, sep='\t')

