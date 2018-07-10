#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import sys
import seaborn as sns
import sys
import numpy as np
from matplotlib.colors import SymLogNorm

def log(message):
    print(message, file=sys.stderr, flush=True)

if __name__ == "__main__":
    data_path = sys.argv[1]
    out_path = sys.argv[2]
    data = pd.read_table(data_path, names=['contig_id', 'position',
                                           'total_depth', 'cosine_similarity'])
    data['log10_total_depth'] = np.log10(data['total_depth'])
    plt.hist2d('log10_total_depth', 'cosine_similarity', data=data,
               bins=1000, cmap='Reds', norm=SymLogNorm(linthresh=1))
    plt.colorbar()
    plt.savefig(out_path)
