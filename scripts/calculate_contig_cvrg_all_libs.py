#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import sys
import numpy as np

if __name__ == "__main__":
    depth = pd.read_table(sys.argv[1],
                          names=['library_id', 'contig_id',
                                 'position', 'depth'])
    length = pd.read_table(sys.argv[2], index_col=['contig_id'], squeeze=True)
    data = (depth.groupby(['contig_id', 'library_id'])
                 .apply(lambda d: d.depth.sum())
                 .unstack('library_id')
                 .apply(lambda x: x / length))
    data.to_csv(sys.stdout, sep='\t')
