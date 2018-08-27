#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import sys

def union_of_intervals(interval_list):
    interval_list = [(min(start, end), max(start, end)) for start, end in interval_list]
    out = []
    for start, end in sorted(interval_list):
        if out and out[-1][1] >= start - 1:
            out[-1][1] = max(out[-1][1], end)
        else:
            out.append([start, end])
    return out

def total_covered_length(interval_list):
    reduced = union_of_intervals(interval_list)
    tally = 0
    for interval in reduced:
        tally += interval[1] - interval[0]
    return tally

def total_covered_in_interval(interval, other_intervals):
    pre_length = total_covered_length(other_intervals)
    post_length = total_covered_length(other_intervals + [interval])
    interval_length = total_covered_length([interval])
    return interval_length - (post_length - pre_length)

def select_domains(df, max_frac_overlap=0.5):
    selected = []
    selected_intervals = []
    for idx, hit in df.sort_values('score', ascending=False).iterrows():
        hit_length = hit.right - hit.left
        overlapping_length = total_covered_in_interval((hit.left, hit.right), selected_intervals)
        if overlapping_length < hit_length * max_frac_overlap:
            selected.append(idx)
            selected_intervals.append((hit.left, hit.right))
    return df.loc[selected]

if __name__ == "__main__":
    domains_path = sys.argv[1]
    max_frac_overlap = float(sys.argv[2])

    data = pd.read_table(domains_path,
                         names=['feature_id', 'domain_id',
                                'score', 'left', 'right'])
    out = (data.groupby('feature_id')
               .apply(lambda d: select_domains(d, max_frac_overlap))
               .reset_index(drop=True)
               .sort_values(['feature_id', 'left']))
    out[['feature_id', 'domain_id', 'left', 'right']].to_csv(sys.stdout, sep='\t', index=False, header=False)
