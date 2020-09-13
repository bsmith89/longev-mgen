#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from Bio.SeqIO import index as index_seqs
import sys
from collections import defaultdict
import argparse


def log(message):
    print(message, file=sys.stderr, flush=True)


def idxwhere(x):
    return x[x].index


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--fasta-output-template", type=str)
    p.add_argument("--tsv-output", type=str)
    p.add_argument("--list-output", type=str)
    p.add_argument("--min-frac", type=float, default=1.0)
    p.add_argument(
        "hits", type=str, nargs="+", metavar="GENOME_ID:TABLE:FASTA"
    )

    args = p.parse_args()

    input_paths = {}
    for arg in args.hits:
        genome_id, hmmer_path, fasta_path = arg.strip().split(":")
        input_paths[genome_id] = (hmmer_path, fasta_path)

    data = []
    for genome_id in input_paths:
        hmmer_path, _ = input_paths[genome_id]
        d = pd.read_table(
            hmmer_path, names=["feature_id", "model_id", "score"]
        )
        d["genome_id"] = genome_id
        data.append(d)
    data = pd.concat(data)

    model_by_genome_counts = (
        data.groupby(["genome_id"])
        .apply(lambda d: d.model_id.value_counts())
        .unstack(fill_value=0)
    )

    is_single_copy = (model_by_genome_counts > 1).sum() == 0
    is_ubiquitous = (model_by_genome_counts == 1).mean() >= args.min_frac

    shared_single_copy_models = idxwhere(is_single_copy & is_ubiquitous)

    if args.list_output:
        with open(args.list_output, "w") as f:
            for model_id in shared_single_copy_models:
                print(model_id, file=f)

    if args.tsv_output:
        (
            data[["genome_id", "model_id", "feature_id"]].to_csv(
                args.tsv_output, sep="\t", header=False, index=False
            )
        )

    if args.fasta_output_template:
        fasta_indexes = {}
        for genome_id in input_paths:
            fasta_path = input_paths[genome_id][1]
            fasta_indexes[genome_id] = index_seqs(fasta_path, "fasta")
        for model_id in shared_single_copy_models:
            with open(args.fasta_output_template.format(model_id), "w") as f:
                for _, x in data[data.model_id == model_id].iterrows():
                    feature_id = x["feature_id"]
                    rec = fasta_indexes[x["genome_id"]][feature_id]
                    print(f">{rec.id}\n{rec.seq}", file=f)
