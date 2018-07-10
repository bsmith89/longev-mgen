#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import sys
import seaborn as sns

if __name__ == "__main__":
    data_path = sys.argv[1]
    out_path = sys.argv[2]
    data = pd.read_table(data_path, names=['contig_id', 'position', 'total_depth', 'cosine_similarity'])
    plt.scatter('total_depth', 'cosine_similarity', data=data, s=0.2, alpha=0.5)
    plt.xscale('log')
    sns.kdeplot(data.total_depth, data.cosine_similarity)
    plt.savefig(out_path)
