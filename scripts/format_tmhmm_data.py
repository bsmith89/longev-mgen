#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import sys

if __name__ == "__main__":
    data = pd.read_table(sys.argv[1], names=['feature_id', '_', 'location', 'start', 'end'], comment='#')
    out = data.groupby(['feature_id', 'location']).apply(len).unstack('location', fill_value=0)
    out['TMhelix'].astype(int).to_csv(sys.stdout, sep='\t', header=False)

