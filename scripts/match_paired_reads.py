#!/usr/bin/env python3
"""Take an unsorted interleaved paired-end FASTQ file and match up reads."""

from Bio.SeqIO import parse, write
import sys


def get_seq_id(line):
    return line.split('\t')[0]

if __name__ == "__main__":
    unmatched_reads = {}
    for line in sys.stdin:
        pair_id, _ = line.split('\t', 1)
        if pair_id in unmatched_reads:
            print(unmatched_reads[pair_id], end='')
            print(line, end='')
            del unmatched_reads[pair_id]
        else:
            unmatched_reads[pair_id] = line
    print('{} unmatched reads'.format(len(unmatched_reads)),
        file=sys.stderr)
    if len(sys.argv) > 1:  # User provided filename for unmatched read IDs
        with open(sys.argv[1], 'w') as handle:
            print('Opening {} to output unmatched reads'.format(sys.argv[1]),
                file=sys.stderr)
            for read in unmatched_reads.values():
                print(read, file=handle, end='')
    else:
        print('No unmatched read filename provided.', file=sys.stderr)


