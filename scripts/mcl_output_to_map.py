#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

if __name__ == "__main__":
    n_clusts = 0
    with open(sys.argv[1]) as handle:
        for line in handle:
            n_clusts += 1
    padding = len(str(n_clusts + 1))
    with open(sys.argv[1]) as handle:
        for i, line in enumerate(handle, start=1):
            orf_id_list = line.strip().split('\t')
            if len(orf_id_list) == 1:
                continue
            for orf_id in line.strip().split('\t'):
                index_str = str(i)
                padding_str = '0' * (padding - len(index_str))
                print(orf_id.strip(), 'Opf' + padding_str + index_str, sep='\t')
