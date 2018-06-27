#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Combine TSVs without headers, potentially including a custom field for each table.

USAGE: combine_tables.py [custom1=]path_to_table_1 [[custom2=]path_to_table_2 [[custom3=]path_to_table_3 [...]]]

"""

import pandas as pd
import sys

if __name__ == "__main__":
    table_args = sys.argv[1:]
    tables = []
    n_cols = None
    for table_str in table_args:
        *custom, path = table_str.split('=')
        this_table = pd.read_table(path, header=None)
        if n_cols is not None:
            assert n_cols == this_table.shape[1], \
                    f"All tables must have the same number of columns {n_cols}."
        else:
            n_cols = this_table.shape[1]
        if custom:
            this_table['custom'] = custom
        tables.append(this_table)
    out = pd.concat(tables)
    if 'custom' in out.columns:  # Make sure that custom is the last column.
        out = out[list(range(n_cols)) + ['custom']]
    out.to_csv(sys.stdout, sep='\t', header=False, index=False)

