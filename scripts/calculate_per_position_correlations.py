#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""USAGE: calculate_per_position_correlations.py library-depth.tsv [library.list] > output.tsv"""

import pandas as pd
import tqdm import tqdm
from scipy.stats import pearsonr

if __name__ == "__main__":
    tqdm.pandas()  # Register tqdm as a pandas progress indicator.
    depth_path = sys.argv[1]
    # Load library list.
    libraries = []
    if len(sys.argv) > 2:
        with open(sys.argv[2]) as library_handle:
            libraries = [line.strip() for line in library_handle]
    else:
        print("No libraries loaded.  All libraries will be considered.".format(depth_path), file=sys.stderr)
    print("Loading depth data from {}.".format(depth_path), file=sys.stderr)
    depth = (pd.read_table(depth_path,
                           names=['library_id', 'contig_id', 'position', 'depth'],
                           dtype=[str, str, int, int],
                           index_col=['contig_id', 'position', 'library_id'])
               .depth.unstack('library_id', fill_value=0))
    # TODO: Filter by library before unstacking.
    if libraries:
        depth = depth[libraries]
    print("Calculating total depth in each library.".format(path), file=sys.stderr)
    library_depth = depth.sum()
    print("Calculating correlation for each position.".format(path), file=sys.stderr)
    correlation = depth.progress_apply(lambda x: pearsonr(library_depth, x)[0])
    correlation.name == 'depth_correlation'
    print("Outputting correlation data.".format(path), file=sys.stderr)
    correlation.dropna().to_csv(sys.stdout, sep='\t')
