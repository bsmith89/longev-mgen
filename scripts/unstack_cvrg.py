#!/usr/bin/env python3

import sys
import pandas as pd

if __name__ == '__main__':
    cvrg = pd.read_table(sys.argv[1],
                         index_col=['library_id', 'contig_id']).coverage
    float_format = sys.argv[2]
    data_t = cvrg.unstack('library_id').fillna(0)
    data_t.to_csv(sys.stdout, sep='\t', float_format=float_format)
