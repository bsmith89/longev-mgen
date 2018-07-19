#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import sys
from Bio.SeqIO import index as index_seqs

if __name__ == "__main__":
    table_path = sys.argv[1]
    seq_paths = dict(arg.strip().split(':') for arg in sys.argv[2:])

    seq_indices = {}
    for gene_id in seq_paths:
        seq_indices[gene_id] = index_seqs(seq_paths[gene_id], 'fasta')

    genome_table = pd.read_table(table_path, names=['genome_id', 'model_id', 'orf_id'])
    gene_list = genome_table.model_id.unique()
    for genome_id, d1 in genome_table.groupby('genome_id'):
        sequence = ''
        for gene_id, d2 in d1.set_index('model_id').loc[gene_list].iterrows():
            sequence += str(seq_indices[gene_id][d2.orf_id].seq)
        print(f'>{genome_id}\n{sequence}')

