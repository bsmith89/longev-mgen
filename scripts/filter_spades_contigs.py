#!/usr/bin/env python3

from Bio.SeqIO import parse, write
import pandas as pd
import numpy as np
import sys

def weighted_median(values, weights=None):
    if weights is None:
        weights = np.ones_like(values)
    d = pd.DataFrame({'value': values, 'weight': weights})
    median_cumweight = d['weight'].sum() / 2
    d = d.sort_values('value')
    d['cumweight'] = d['weight'].cumsum()
    return d[d['cumweight'] > median_cumweight].head(1)['value'].values[0]

if __name__ == "__main__":
    length_thresh = int(sys.argv[1])
    cov_thresh = float(sys.argv[2])  # Fraction of median coverage

    recs = {rec.id: rec for rec in parse(sys.stdin, 'fasta')}
    tabular = []
    for rec_id in recs:
        _, i, _, length, _, cov = rec_id.split('_')
        length = int(length)
        cov = float(cov)
        tabular.append((rec_id, length, cov))
    data = pd.DataFrame(tabular, columns=['rec_id', 'length', 'coverage']).set_index('rec_id')
    print("Processed {} sequences".format(data.shape[0]), file=sys.stderr)

    data_len_filt = data[lambda x: x.length > length_thresh]
    print("{} sequences pass the length threshold ({})".format(data_len_filt.shape[0], length_thresh), file=sys.stderr)
    median_cov = weighted_median(data_len_filt.coverage, data_len_filt.length)
    print("Median coverage of {}".format(median_cov), file=sys.stderr)
    data['coverage_norm'] = data.coverage / median_cov
    data_filt = data[lambda x: (x.coverage_norm > cov_thresh) &
                               (x.length > length_thresh)]
    print("{} sequences pass the coverage threshold ({})".format(data_filt.shape[0], cov_thresh), file=sys.stderr)

    for rec_id in data_filt.index:
        write(recs[rec_id], sys.stdout, 'fasta')



