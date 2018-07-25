#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import sys
from sklearn.cluster import AgglomerativeClustering


def matrix_from_sparse(data, diag=0, offdiag=1):
    """Construct a (dis)similarity matrix from a sparse representation.

    Default arguments are for a dissimilarity index that is bound between 0
    and 1.

    """
    # Un-sparsify
    dmatrix = data.unstack()
    # Fill out the axes
    orfs = sorted(list(set(dmatrix.index) | set(dmatrix.columns)))
    dmatrix = dmatrix.reindex(index=orfs, columns=orfs)
    # Mirror the data over the diagonal
    dmatrix = dmatrix.fillna(dmatrix.T)
    # Set the diagonal
    np.fill_diagonal(dmatrix.values, 0)
    # Fill any other gaps
    dmatrix = dmatrix.fillna(1)
    return dmatrix


if __name__ == "__main__":
    data = pd.read_table(sys.argv[1],
                     names=['qseqid', 'sseqid', 'dissimilarity'], index_col=['qseqid', 'sseqid'])
    dmatrix = matrix_from_sparse(data.dissimilarity)

    ac = AgglomerativeClustering(n_clusters=int(sys.argv[2]),
                                 affinity='precomputed',
                                 linkage='average').fit(dmatrix)
    cluster = pd.DataFrame({'cluster': ac.labels_}, index=dmatrix.index)

    cluster_size = cluster.reset_index().groupby('cluster').qseqid.count()
    is_singleton_cluster = (cluster_size == 1)
    is_true_cluster = (cluster_size > 1)
    true_clusters = cluster_size[is_true_cluster].index
    in_true_clusters = cluster.cluster.isin(true_clusters)
    n_clusters = len(true_clusters)
    n_in_clusters = in_true_clusters.sum()
    print(f'Found {n_clusters} clusters with {n_in_clusters} total CDS.', file=sys.stderr)

    largest_cluster = cluster_size.idxmax()
    largest_members = cluster[lambda x: x.cluster == largest_cluster].index
    largest_n_members = len(largest_members)
    largest_dmatrix = dmatrix.loc[largest_members, largest_members]
    largest_max_dist = largest_dmatrix.max().max()
    print(f'Largest cluster has {largest_n_members} members with a maximum dissimilarity of {largest_max_dist}.', file=sys.stderr)

    out = cluster.loc[lambda x: x.cluster.isin(true_clusters)]  # Filter out singletons
    n_digits = int(np.ceil(np.log10(int(sys.argv[2]) + 1)))
    template = 'Opu{:0' + str(n_digits) + '}'
    out.cluster = out.cluster.apply(lambda n: template.format(n))
    out.to_csv(sys.stdout, sep='\t', header=False)
