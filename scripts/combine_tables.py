#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Combine TSVs without headers, potentially including a custom field for each table.

USAGE: combine_tables.py [custom1=]path_to_table_1 [[custom2=]path_to_table_2 [[custom3=]path_to_table_3 [...]]]

"""

import pandas as pd
import sys

if __name__ == "__main__":
    table_args = sys.argv[1:]
    n_cols = None
    for table_str in table_args:
        *custom, path = table_str.split('=')
        table = pd.read_table(path, header=None)

        # Defensively check column count.
        if n_cols is not None:
            assert n_cols == table.shape[1], \
                    f"All tables must have the same number of columns {n_cols}."
        else:
            n_cols = table.shape[1]

        # Add a custom value if it exists.
        if custom:
            table['custom'] = custom[0]

        # Defensively assign column order.
        if 'custom' in table.columns:
            table = table[list(range(n_cols)) + ['custom']]
        else:
            table = table[list(range(n_cols))]

        table.to_csv(sys.stdout, sep='\t', header=False, index=False)

