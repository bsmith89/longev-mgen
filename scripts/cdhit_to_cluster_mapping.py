#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys


def chunk_lines_by_cluster(lines):
    clustnum=None
    outlines = []
    for line in lines:
        if line.startswith('>Cluster '):
            if outlines:
                yield clustnum, outlines
            outlines = []
            clustnum = int(line.strip().split()[1])
        else:
            outlines.append(line)
    else:
        yield clustnum, outlines


def cluster_lines_to_details(lines):
    seqlist = []
    leadname = None
    for line in lines:
        seqnum, length, seqidstr, *rest = line.strip().split()
        seqid = seqidstr[1:-3]
        if rest[0] == '*':
            leadname = seqid
        seqlist.append(seqid)
    assert leadname is not None
    return leadname, seqlist


if __name__ == "__main__":
    with open(sys.argv[1]) as f:
        for clustnum, chunk in chunk_lines_by_cluster(f):
            leadname, seqids = cluster_lines_to_details(chunk)
            for seq in seqids:
                print(seq, leadname, sep='\t')
