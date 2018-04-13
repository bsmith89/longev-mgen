#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio.SearchIO import parse
import sys

if __name__ == "__main__":
    print('hit_id', 'hit_description', 'query_id', 'hit_start', 'hit_end',
          'query_start', 'query_end', 'bitscore', 'included', sep='\t')
    recs = parse(sys.stdin, 'hmmer3-text')
    for model_hits in recs:
        for hsp in model_hits.hsps:
            print(hsp.hit_id, hsp.hit_description,
                  hsp.query_id,
                  hsp.hit_start, hsp.hit_end,
                  hsp.query_start, hsp.query_end,
                  hsp.bitscore, int(hsp.is_included),
                  sep='\t')
