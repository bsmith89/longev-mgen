#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import sys


def mean_matrix(arr):
    "Construct N-by-N matrix where Y_ij  is the mean of X_i and X_j."
    return np.apply_along_axis(lambda x: (x + arr) / 2,
                               axis=1,
                               arr=np.reshape(arr, (-1, 1)))


if __name__ == "__main__":
    data = pd.read_table(sys.argv[1],
                         names=['qseqid', 'sseqid', 'pident', 'length',
                                 'mismatch','gapopen', 'qstart', 'qend',
                                 'sstart', 'send', 'evalue', 'bitscore'])
    similarity = data.groupby(['qseqid', 'sseqid']).bitscore.sum()
    smatrix = similarity.unstack()  # Make into a matrix
    smatrix = (smatrix + smatrix.T) / 2  # Make symmetrical
    assert (smatrix.columns == smatrix.index).all()
    smatrix = smatrix / mean_matrix(np.diag(smatrix))  # Normalize to mean
    dmatrix = 1 - smatrix  # Convert to dissimilarity
    (dmatrix.stack(dropna=True).reset_index()
            [lambda x: (x.qseqid > x.sseqid)]
            .to_csv(sys.stdout, sep='\t', header=False, index=False)
    )
