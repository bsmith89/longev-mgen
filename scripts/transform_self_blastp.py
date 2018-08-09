#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import sys

data = pd.read_table(sys.argv[1],
                     names=['qseqid', 'sseqid', 'pident', 'length',
                            'mismatch','gapopen', 'qstart', 'qend',
                            'sstart', 'send', 'evalue', 'bitscore'])
similarity = data.groupby(['qseqid', 'sseqid']).bitscore.sum().reset_index()
similarity = pd.merge(similarity, similarity, how='outer',
                      left_on=['qseqid', 'sseqid'], right_on=['sseqid', 'qseqid'])
similarity[['bitscore_x', 'bitscore_y']].fillna(0, inplace=True)
similarity.set_index(['qseqid_x', 'sseqid_x'], inplace=True)
similarity['mean_bitscore'] = similarity.apply(lambda x: (x.bitscore_x + x.bitscore_y) / 2, axis='columns')
similarity['mean_bitscore'].reset_index()[lambda x: x.qseqid_x < x.sseqid_x].to_csv(sys.stdout, index=False, sep='\t')
