#!/usr/bin/env python3
"""python3 plot_contig_ordination.py features_file.tsv bins_file.tsv [bin1 [bin2 [...]]] output_name.pdf"""

import pandas as pd
import matplotlib.pyplot as plt
import sys
from Bio.SeqIO import parse
from sklearn.decomposition import PCA
import logging

logger = logging.getLogger(__file__)

if __name__ == "__main__":
    logging.basicConfig(
            level=logging.DEBUG,
            format='%(message)s')

    tetra = pd.read_table(sys.argv[1], index_col='contig_id')
    tetra = tetra[[str(i) for i in range(136)]]
    logger.debug('Feature data loaded.')
    pca = PCA(0.9, copy=False, random_state=1).fit_transform(tetra)
    pca = pd.DataFrame(pca,
                       columns=['x{:02d}'.format(int(i))
                                for i in range(pca.shape[1])],
                       index=tetra.index)
    logger.debug('PCA calculated.')

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(pca.x00, pca.x01, s=1, alpha=0.2, color='gray', rasterized=True,
               label='_nolegend_')
    logger.debug('Background points plotted.')

    # TODO: Correct the column names for .bins.tsv files.
    bins = pd.read_table(sys.argv[2], index_col='contig_id',
                         skiprows=1, dtype={'bin_id': str},
                         names=['contig_id', 'bin_id'])
    logger.debug('Bin mapping loaded.')
    for bin_id in sys.argv[3:-1]:
        contig_names = list(bins[bins.bin_id == bin_id].index)
        pca_bin = pca.loc[contig_names]
        ax.scatter(pca_bin.x00, pca_bin.x01, s=4, lw=1, label=bin_id, alpha=0.5)
    ax.legend()
    ax.set_xticks([])
    ax.set_yticks([])
    fig.savefig(sys.argv[-1], dpi=240)
    logger.debug('Figure saved to {}'.format(sys.argv[-1]))
