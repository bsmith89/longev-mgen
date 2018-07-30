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
               bins=[np.linspace(0, data.log10_total_depth.max(), 1000),
                     np.linspace(0.5, 1, 1000)],
               cmap='Reds', norm=SymLogNorm(linthresh=1))
    cb = plt.colorbar()
    plt.xlabel("Log10(depth)")
    plt.ylabel("CosSim(depth, trusted_depth)")
    cb.set_label("Count")
    plt.savefig(out_path)
