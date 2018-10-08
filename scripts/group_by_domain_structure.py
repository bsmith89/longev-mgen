#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import sys


if __name__ == "__main__":
    data = pd.read_table(sys.argv[1],
                         names=['feature_id', 'domain_id', 'score', 'left', 'right'])
    (data.sort_values(['feature_id', 'left'])
         .groupby('feature_id').domain_id.apply(lambda x: ':'.join(x))
         .reset_index()
         .to_csv(sys.stdout, sep='\t', header=False, index=False))
