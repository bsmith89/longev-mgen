#!/usr/bin/env python3
"""
Cluster contigs using sklearn.

"""

import sys
import argparse
import logging
from contextlib import redirect_stdout

import pandas as pd
import numpy as np
from sklearn.mixture import GaussianMixture

logger = logging.getLogger(__name__)


def load_data(data_path, length_path, min_length):
    # Load from file.
    data = pd.read_csv(args.data_path, index_col='contig_id')
    length = pd.read_table(args.length_path, index_col='contig_id').length
    # Align indices.
    length = length.loc[data.index]
    # Length filter.
    data = data[length >= min_length]
    length = length.loc[data.index]
    logger.info('Loaded {} contigs with {} bp total.'
                    .format(len(length), length.sum())
               )
    return data, length

def cluster_optics(data, threads=1):
    from sklearn.cluster import OPTICS
    logger.info('Fitting an OPTICS model.')
    model = OPTICS(n_jobs=args.threads)
    model.fit(data)
    logger.info('Finished clustering.')
    logger.info('Assigning bins to contigs.')
    cluster = pd.Series(model.labels_, index=data.index, name='cluster')
    return cluster


def _fit_gmm(data, max_nbins):
    model = GaussianMixture(max_nbins, covariance_type='full',
                            verbose=2, verbose_interval=1)
    with redirect_stdout(sys.stderr):
        model.fit(data)
    if model.converged_:
        logger.info('Finished clustering.')
    else:
        logger.warning('Convergence not achieved.  Results may not be correct.')
    return model


def _label_gmm(model, data, pthresh=0.5):
    logger.info(f'Calculating posterior probabilities.')
    probs = model.predict_proba(data)
    logger.info(f'Assigning bins to contigs where p > {pthresh}.')
    best_guess = np.argmax(probs, axis=1)
    # Are any of the probs > pthresh?
    mask = (probs > pthresh).sum(axis=1).astype(bool)
    cluster = pd.Series(best_guess.astype(int), index=data.index, name='cluster')
    return cluster[mask]


def cluster_gmm(data, max_nbins):
    logger.info('Fitting a GMM model using EM.')
    logger.info(f'max_nbins={max_nbins}')
    model = _fit_gmm(data, max_nbins)
    cluster = _label_gmm(model, data)
    return cluster


def cluster_gmmss(data, length, max_nbins, frac):
    logger.info('Fitting a GMM model with a subsample of data using EM.')
    logger.info(f'frac={frac}, max_nbins={max_nbins}')
    model = _fit_gmm(data.sample(frac=frac, weights=length), max_nbins)
    cluster = _label_gmm(model, data)
    return cluster


def filt_by_total_size(cluster, length, min_bin_size):
    valid_clusters = (length.groupby(cluster).sum() >= min_bin_size)[lambda x: x].index
    filt_cluster = cluster[cluster.isin(valid_clusters)]
    total_clusters = len(valid_clusters)
    logger.info(f'{total_clusters} clusters with more than {min_bin_size} bp')
    return filt_cluster


def summarize_clusters(cluster, length):
    cluster_summary = (length.groupby(cluster).agg(['count', 'sum'])
                             .sort_values('sum', ascending=False))
    cluster_summary.rename(columns={'count': 'n_contigs',
                                    'sum': 'total_length'},
                           inplace=True)
    cluster_summary.rename(index=int, inplace=True)
    logger.info('{} bp in {} valid clusters'
                    .format(cluster_summary.total_length.sum(),
                            cluster_summary.shape[0])
               )
    logger.info('Top 10 clusters:\n\n{}'.format(cluster_summary.head(10)))
    return cluster_summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('data_path', metavar='CONCOCT_PCA_TRANSFORM_CSV')
    parser.add_argument('length_path', metavar='CONTIG_LENGTH_TSV')
    parser.add_argument('--min-length', type=int, default=0)
    parser.add_argument('--min-bin-size', type=int, default=100000,
                        help=('the minimum size (in bp) of a cluster'
                              ' (default: %(default)s)'))
    parser.add_argument('--verbosity', default='DEBUG')
    parser.add_argument('--summary',
                        help='Print summary of clusters to path')

    # Add subparsers for different methods.
    subparsers = parser.add_subparsers(title='clustering methods', dest='method')
    gmm_p = subparsers.add_parser('gmm')
    gmmss_p = subparsers.add_parser('gmm-ss', help='GMM on subsample of data')
    optics_p = subparsers.add_parser('optics')

    # Add arguments to subparsers.
    optics_p.add_argument('--threads', type=int, default=1)
    for subparser in [gmm_p, gmmss_p]:
        subparser.add_argument('--max-nbins', type=int, default=1000)
        subparser.add_argument('--prob-min', type=float, default=0.0,
                               help='posterior probability minimum to be assigned to a bin')
    gmmss_p.add_argument('--frac', type=float, default=0.1, help='fraction of data used to fit initial GMM')

    args = parser.parse_args()

    logging.basicConfig(format='%(message)s',
                        level=getattr(logging, args.verbosity.upper()))

    data, length = load_data(args.data_path, args.length_path, args.min_length)

    if args.method == 'optics':
        all_clusters = cluster_optics(data, threads=args.threads)
    elif args.method == 'gmm':
        all_clusters = cluster_gmm(data, max_nbins=args.max_nbins)
    elif args.method == 'gmm-ss':
        all_clusters = cluster_gmmss(data, length,
                                     max_nbins=args.max_nbins, frac=args.frac)
    else:
        raise NotImplementedError("The {args.method} clustering method has not yet been implemented.")

    filt_clusters = filt_by_total_size(all_clusters, length,
                                       min_bin_size=args.min_bin_size)
    cluster_summary = summarize_clusters(filt_clusters, length)

    if args.summary:
        cluster_summary.to_csv(args.summary, sep='\t')
    filt_clusters.to_csv(sys.stdout, sep='\t', header=True)

