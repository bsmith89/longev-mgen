#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from Bio.SeqIO import index as index_seqs
import pandas as pd
import numpy as np

if __name__ == "__main__":
    table_path = sys.argv[1]
    align_path = dict(s.split(":") for s in sys.argv[2:])
    model_list = list(align_path.keys())

    gene_table = pd.read_table(
        table_path, names=["genome_id", "model_id", "feature_id"]
    )
    genome_list = gene_table.genome_id.unique()

    align_index = {}
    for model_id in model_list:
        align_index[model_id] = index_seqs(align_path[model_id], "fasta")

    align_len = {}
    for model_id, d in gene_table[
        gene_table.model_id.isin(model_list)
    ].groupby("model_id"):
        align_len[model_id] = len(
            align_index[model_id][d.iloc[0].feature_id].seq
        )
        print(
            model_id,
            align_path[model_id],
            align_len[model_id],
            file=sys.stderr,
        )

    for genome_id, d in gene_table.groupby("genome_id"):
        d = (
            d[d.model_id.isin(model_list)]
            .set_index("model_id")
            .reindex(model_list)
        )
        d.feature_id.replace({float('nan'): None}, inplace=True)

        sequence = ""
        length = None
        for model_id, x in d.iterrows():
            if x.feature_id is None:
                seq = "-" * align_len[model_id]
            else:
                seq = align_index[model_id][x.feature_id].seq
                assert len(seq) == align_len[model_id]
            sequence += seq
        if length is not None:
            assert len(sequence) == length
        print(genome_id, len(sequence), file=sys.stderr)
        print(f">{genome_id}\n{sequence}")
