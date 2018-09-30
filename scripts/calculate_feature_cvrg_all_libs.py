#!/usr/bin/env python3

import pandas as pd
import sys
import numpy as np

def mean_feature_depth(depth, feature):
    """Mean depth of each feature."""
    length = feature.left - feature.right + 1
    feature_intervals = pd.IntervalIndex.from_arrays(feature.left, feature.right)
    sum_depth = (depth.groupby(['sequence_id',
                                pd.cut(depth.position,
                                       feature_intervals)])
                      .depth.sum()
                      .reindex(feature.index)
                )
    return sum_depth.divide(length)

if __name__ == "__main__":
    depth = pd.read_table(sys.argv[1],
                          names=['library_id', 'sequence_id',
                                 'position', 'depth'])
    feature = pd.read_table(sys.argv[2],
                            names=['feature_id', 'sequence_id',
                                   'left', 'right'])
    feature_sum_depth = (depth.merge(feature, on=['sequence_id'])
                              [lambda x: x.position.between(x.left, x.right)]
                              .groupby(['library_id', 'feature_id'])
                              .depth.sum()
                        )
    feature_length = feature.set_index('feature_id').apply(lambda x: np.abs(x.right - x.left) + 1, axis='columns')
    feature_mean_depth = (feature_sum_depth
                              .groupby(level='library_id')
                              .apply(lambda x: x / feature_length)
                         )
    feature_mean_depth.to_csv(sys.stdout, sep='\t')
    # (depth.groupby('library_id')
    #       .apply(lambda x: mean_feature_depth(x, feature))
    #       .to_csv(sys.stdout, sep='\t')
    # )
