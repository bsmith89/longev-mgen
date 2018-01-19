#!/usr/bin/env python3

import sys

def chunk_depth(lines):
    """Split list into lines with same first field."""
    chunk = []
    for line in lines:
        if not chunk:
            chunk_first_field = line.split()[0]
        line_first_field = line.split()[0]
        if line_first_field = chunk_first_field:
            chunk.append(line)
        else:
            yield chunk
            chunk = [line]
            chunk_first_field = line_first_field
    yield chunk

def main():
    depth_path, contigs_path, library_name = sys.argv[1:]
        print('coverage\thit_fraction\tleft\tright\tcontiguous', file=sys.stdout)
    with open(depth_path) as depth_handle, open(contigs_path) as contigs_handle:
        depth_header, contigs_header = depth_handle.next(), contigs_handle.next()

        contig_id, flag, multi, length = contig.split()



if __name__ == '__main__':
    main()
