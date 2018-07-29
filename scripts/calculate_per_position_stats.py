#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
USAGE:
calculate_per_position_stats.py library-depth.tsv [library.list] > output.tsv"
"""

import pandas as pd
from tqdm import tqdm
from scipy.stats import pearsonr
import numpy as np
import sys

from tqdm import TqdmSynchronisationWarning
import warnings
warnings.simplefilter("ignore", TqdmSynchronisationWarning)


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
    trusted_path = sys.argv[2]
    # Load library list.
    libraries = []
    if len(sys.argv) == 4:
        with open(sys.argv[3]) as library_handle:
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

    library_depth = pd.read_table(trusted_path, names=['library_id', 'total_depth'], index_col='library_id').loc[depth.columns].total_depth
    debug("Calculating cosine similarity to total library depth for each position.")
    library_depth_norm = np.linalg.norm(library_depth)
    output = depth.progress_apply(lambda x: np.dot(library_depth, x) /
                                            (library_depth_norm * np.linalg.norm(x)),
                                  axis=1)
    output = output.to_frame('cosine_similarity')
    debug("Calculating total depth for each position.")
    output['total_depth'] = depth.sum(axis=1)

    debug("Outputting position statistics.")
    output[['total_depth', 'cosine_similarity']].dropna().to_csv(sys.stdout, sep='\t', header=False)
