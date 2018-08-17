#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

if __name__ == "__main__":
    with open(sys.argv[1]) as handle:
        for line in handle:
            if line.startswith('>Feature '):
                sequence_id = line.strip()[9:]
            elif not line.startswith('\t'):
                left, right, *_ = line.strip().split('\t')
                region_has_previous_locus_tag = False
            elif line.strip().startswith('locus_tag'):
                assert not region_has_previous_locus_tag
                feature_id = line.strip().split('\t')[1]
                region_has_previous_locus_tag = True
                print(feature_id, sequence_id, left, right, sep='\t')
