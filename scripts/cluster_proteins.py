#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import sys
from sklearn.cluster import AgglomerativeClustering


def binary_operator_matrix(arr, func):
    """Construct N-by-N matrix where Y_ij is the result of func(X_i and X_j).

    *func(arr, axis)* : a function like numpy.mean, that takes a keyword *axis*.
        np.max, np.min, etc. all work.

    """
    return func(np.stack([np.repeat(arr, len(arr)).reshape(-1, len(arr), order='F'),
                          np.repeat(arr, len(arr)).reshape(-1, len(arr), order='C')]),
                axis=0)


if __name__ == "__main__":

    input_path = sys.argv[1]
    n_clusters = int(sys.argv[2])

    data = pd.read_table(input_path,
                         names=['qseqid', 'sseqid', 'pident', 'length',
                                'mismatch','gapopen', 'qstart', 'qend',
                                'sstart', 'send', 'evalue', 'bitscore'])
    similarity = data.groupby(['qseqid', 'sseqid']).bitscore.sum()
    smatrix = similarity.unstack()  # Make into a matrix
    smatrix = (smatrix + smatrix.T) / 2  # Make symmetrical
    assert (smatrix.columns == smatrix.index).all()
    smatrix = smatrix / binary_operator_matrix(np.diag(smatrix), func=np.max)  # Normalize to max
    dmatrix = 1 - smatrix  # Convert to dissimilarity
    del smatrix
    dmatrix = dmatrix.fillna(dmatrix.T)
    # Set the diagonal
    np.fill_diagonal(dmatrix.values, 0)
    # Fill any other gaps
    dmatrix = dmatrix.fillna(1)

    ac = AgglomerativeClustering(n_clusters=n_clusters,
                                 affinity='precomputed',
                                 linkage='average').fit(dmatrix)
    cluster = pd.DataFrame({'cluster': ac.labels_}, index=dmatrix.index)

    cluster_size = cluster.reset_index().groupby('cluster').qseqid.count()
    is_singleton_cluster = (cluster_size == 1)
    is_true_cluster = (cluster_size > 1)
    true_clusters = cluster_size[is_true_cluster].index
    in_true_clusters = cluster.cluster.isin(true_clusters)
    count_clusters = len(true_clusters)
    n_in_clusters = in_true_clusters.sum()
    print(f'Found {count_clusters} clusters with {n_in_clusters} total CDS.', file=sys.stderr)

    largest_cluster = cluster_size.idxmax()
    largest_members = cluster[lambda x: x.cluster == largest_cluster].index
    largest_n_members = len(largest_members)
    largest_dmatrix = dmatrix.loc[largest_members, largest_members]
    largest_max_dist = largest_dmatrix.max().max()
    print(f'Largest cluster has {largest_n_members} members with a maximum dissimilarity of {largest_max_dist}.', file=sys.stderr)

    out = cluster.loc[lambda x: x.cluster.isin(true_clusters)]  # Filter out singletons
    n_digits = int(np.ceil(np.log10(n_clusters + 1)))
    template = 'Opu{:0' + str(n_digits) + '}'
    out.cluster = out.cluster.apply(lambda n: template.format(n))
    out.to_csv(sys.stdout, sep='\t', header=False)
