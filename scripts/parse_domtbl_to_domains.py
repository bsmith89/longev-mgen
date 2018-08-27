#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import Bio.SearchIO
import pandas as pd
import sys

if __name__ == "__main__":
    domtbl_path = sys.argv[1]

    recs = []
    for query in Bio.SearchIO.parse(domtbl_path, 'hmmsearch3-domtab'):
        for i in range(len(query)):
            hit = query[i]
            for j in range(len(hit)):
                hsp = hit[j]
                for k in range(len(hsp)):
                    frag = hsp[k]
                    recs.append((frag.hit_id, frag.query_id, hsp.bitscore, frag.hit_start, frag.hit_end))
    data = (pd.DataFrame(recs, columns=['feature_id', 'domain_id', 'score',
                                        'left', 'right'])
              .sort_values('feature_id'))
    data[['feature_id', 'domain_id', 'score', 'left', 'right']].to_csv(sys.stdout, sep='\t', index=False, header=False)
