#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
USAGE:
calculate_per_position_correlations.py library-depth.tsv [library.list] > output.tsv"
"""

import pandas as pd
from tqdm import tqdm
from scipy.stats import pearsonr
import sys

# import numpy as np
# from multiprocessing import cpu_count, Parallel
#
# def parallelize(data, func, partitions=None):
#     if partitions is None:
#         partitions = cpu_count()
#     data_split = np.array_split(data, partitions)
#     pool = Pool(cores)
#     data = pd.concat(pool.map(func, data_split))
#     pool.close()
#     pool.join()
#     return data

def debug(s):
    print(s, file=sys.stderr, flush=True)

if __name__ == "__main__":
    tqdm.pandas()  # Register tqdm as a pandas progress indicator.
    depth_path = sys.argv[1]
    # Load library list.
    libraries = []
    if len(sys.argv) == 3:
        with open(sys.argv[2]) as library_handle:
            libraries = [line.strip() for line in library_handle]
    else:
        debug("Library list empty.  All libraries will be considered.")

    debug("Loading depth data from {}.".format(depth_path))
    depth = pd.read_table(depth_path,
                          names=['library_id', 'contig_id', 'position', 'depth'],
                          dtype={'library_id': str, 'contig_id': str,
                                 'position': int, 'depth': int})
    if libraries:
        depth = depth[depth.library_id.isin(libraries)]
    depth.set_index(['library_id', 'contig_id', 'position'], inplace=True)
    debug("Unstacking depth data.")
    depth = depth.depth.unstack('library_id', fill_value=0)

    debug("Calculating total depth in each library.")
    library_depth = depth.sum()
    debug("Calculating correlation for each position.")
    correlation = depth.progress_apply(lambda x: pearsonr(library_depth, x)[0], axis=1)
    correlation = correlation.to_frame('depth_correlation')

    debug("Outputting correlation data.")
    correlation.dropna().to_csv(sys.stdout, sep='\t', header=False)
