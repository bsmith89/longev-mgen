#!/usr/bin/env python3

from Bio.SeqIO import parse, write
import pandas as pd
import sys

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
    median_cov = data_len_filt.coverage.median()
    print("Median coverage of {}".format(median_cov), file=sys.stderr)
    data['coverage_norm'] = data.coverage / median_cov
    data_filt = data[lambda x: (x.coverage_norm > cov_thresh) &
                               (x.length > length_thresh)]
    print("{} sequences pass the coverage threshold ({})".format(data_filt.shape[0], cov_thresh), file=sys.stderr)

    for rec_id in data_filt.index:
        write(recs[rec_id], sys.stdout, 'fasta')



