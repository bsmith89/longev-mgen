#!/usr/bin/env python3
"""Take an unsorted interleaved paired-end FASTQ file and match up reads."""

from Bio.SeqIO import parse, write
import sys


def chunk_fastq(lines):
    chunk = []
    for i, line in enumerate(lines):
        if i % 4 == 0:
            if chunk:
                yield chunk
            chunk = [line]
        else:
            chunk.append(line)
    yield chunk


def get_seq_id(chunk):
    return chunk[0][1:-1]

if __name__ == "__main__":
    unmatched_reads = {}
    for chunk in chunk_fastq(sys.stdin):
        seq_id = get_seq_id(chunk)
        pair_id = seq_id[:-2]
        if pair_id in unmatched_reads:
            read_a = unmatched_reads[pair_id]
            read_b = chunk
            if get_seq_id(read_a) < get_seq_id(read_b):
                print(''.join(read_a), end='')
                print(''.join(read_b), end='')
            else:
                print(''.join(read_b), end='')
                print(''.join(read_a), end='')
            del unmatched_reads[pair_id]
        else:
            unmatched_reads[pair_id] = chunk
    if unmatched_reads:
        print('{} unmatched reads appended to end of output'.format(len(unmatched_reads)),
              file=sys.stderr)
        for chunk in unmatched_reads.values:
            print(''.join(chunk), end='')

