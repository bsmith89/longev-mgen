#!/usr/bin/env python3

import pandas as pd
import sys

def load_list(path):
    with open(path) as handle:
        out = [line.strip() for line in handle]
    return out


def load_cvrg(path):
    out = (pd.read_table(path,
                          names=['library_id', 'feature_id', 'coverage'],
                          index_col=['library_id', 'feature_id'],
                          squeeze=True)
             .unstack('feature_id', fill_value=0)
             .apply(lambda x: x / x.median(), axis=1)
          )
    return out


if __name__ == "__main__":
    cvrg = load_cvrg(sys.argv[1])
    lib_1 = load_list(sys.argv[2])
    lib_2 = load_list(sys.argv[3])
    assert not set(lib_1) & set(lib_2), "Library lists overlap."
    lib_1 = set(lib_1) & set(cvrg.index)
    lib_2 = set(lib_2) & set(cvrg.index)

    cvrg_1 = cvrg.loc[lib_1].mean()
    cvrg_2 = cvrg.loc[lib_2].mean()
    ratio = cvrg_2 / cvrg_1
    ratio.to_csv(sys.stdout, sep='\t', header=False)
