#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import sys
from pandarallel import pandarallel

if __name__ == "__main__":
    pandarallel.initialize(nb_workers=sys.argv[1])

    print("Reading blastp results in", sys.argv[2], file=sys.stderr)
    blastp = pd.read_table(
        sys.argv[2],
        names=[
            "seqA",
            "seqB",
            "_3",
            "_4",
            "_5",
            "_6",
            "_7",
            "_8",
            "_9",
            "_10",
            "_11",
            "score",
        ],
        usecols=["seqA", "seqB", "score"],
    )
    print(blastp.shape[0], "records found", file=sys.stderr)
    blastp[["seqA", "seqB"]] = blastp[["seqA", "seqB"]].where(
        blastp.seqA <= blastp.seqB, blastp[["seqB", "seqA"]]
    )
    # TODO: Consider summing scores here, although that might be a bad idea in
    # cases where the hits are multiplied by repetitive domains in one or both
    # sequences.
    # Can this be corrected by selecting the worse of the two summed scores?
    blastp = blastp.sort_values("score").drop_duplicates(
        ["seqA", "seqB"], keep="last"
    )
    self_score = blastp[lambda x: x.seqA == x.seqB]
    print(self_score.shape[0], "self-hits found", sys.argv[2], file=sys.stderr)
    blastp = blastp.join(
        self_score.set_index("seqA").score, on="seqA", rsuffix="_selfA"
    ).join(self_score.set_index("seqB").score, on="seqB", rsuffix="_selfB")
    blastp["max_self_score"] = blastp[["score_selfA", "score_selfB"]].max(
        axis="columns"
    )
    print(
        "Normalizing scores for {} records in".format(blastp.shape[0]),
        sys.argv[2],
        file=sys.stderr,
    )
    blastp["normalized_score"] = blastp.parallel_apply(
        lambda x: x.score / x.max_self_score, axis="columns"
    )
    print(
        "Outputting",
        blastp.dropna(subset=["normalized_score"]).shape[0],
        "records from",
        sys.argv[2],
        file=sys.stderr,
    )

    (
        blastp[["seqA", "seqB", "normalized_score"]].to_csv(
            sys.stdout, sep="\t", header=False, index=False
        )
    )
