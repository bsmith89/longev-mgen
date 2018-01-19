#!/usr/bin/env python3

import sys

def chunk_depth(handle):
    """Split depth file lines into tables by the contig_id."""
    chunk = []
    for line in handle:
        contig_id, position, depth = line.split()
        if not chunk:
            chunk_contig_id = contig_id
        line_contig_id = contig_id
        if line_contig_id == chunk_contig_id:
            chunk.append((contig_id, int(position), int(depth)))
        else:
            yield chunk
            chunk = [(contig_id, int(position), int(depth))]
            chunk_contig_id = line_contig_id
    yield chunk



def main():
    depth_path, contigs_path, library_id = sys.argv[1:]
    print('library_name', 'contig_id',
          'coverage', 'hit_fraction',
          'left', 'right', 'contiguous',
          sep='\t', file=sys.stdout)

    with open(depth_path) as depth_handle, open(contigs_path) as contigs_handle:
        # Skip the header lines
        depth_header = next(depth_handle)
        contigs_header = next(contigs_handle)
        for depth_table in chunk_depth(depth_handle):
            contig_id = depth_table[0][0]
            contig_id_b = None
            # Scan for the current contig in the metadata
            # (This depends on contigs being in the same order in both files.)
            while contig_id_b != contig_id:
                contig_id_b, flag, multi, length = next(contigs_handle).split()
            length = int(length)
            # Transpose the table, and extract a list of positions and depths.
            _, positions, depths = zip(*depth_table)
            coverage = sum(depths) / length
            num_positions = len(positions)
            hit_fraction = num_positions / length
            left = positions[0]
            right = positions[-1]
            contiguous = int((right - left) + 1 == num_positions)
            print(library_id, contig_id,
                  coverage, hit_fraction,
                  left, right, contiguous,
                  sep='\t', file=sys.stdout)





if __name__ == '__main__':
    main()
