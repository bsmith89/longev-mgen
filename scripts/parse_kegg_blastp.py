#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

if __name__ == "__main__":
    for line in sys.stdin:
        qseqid, sseqid = line.strip().split('\t')
        _, ko_str = sseqid.split('|')
        for ko in ko_str.split(','):
            print(qseqid, ko, sep='\t')

