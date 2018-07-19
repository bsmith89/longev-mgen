#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import sys


if __name__ == "__main__":
    data = pd.read_table(sys.argv[1],
                         names=['domain_id', 'domain_',
                                'query_id', 'query_length',
                                'evalue', 'domain_start',
                                'domain_end', 'query_start',
                                'query_end', 'score'])
    (data.sort_values(['query_id', 'query_start'])
         .groupby('query_id').domain_id.apply(lambda x: ':'.join(x))
         .reset_index()
         .to_csv(sys.stdout, sep='\t', header=False, index=False))
