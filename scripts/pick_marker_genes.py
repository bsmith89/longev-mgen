#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from Bio.SeqIO import index as index_seqs
import sys
from collections import defaultdict
import argparse

def log(message):
    print(message, file=sys.stderr, flush=True)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('--fasta-output-template', type=str, default='marker_genes.{}.fa')
    p.add_argument('--tsv-output', type=str, default='marker_genes.tsv')
    p.add_argument('--list-output', type=str, default='marker_genes.list')
    p.add_argument('hits', type=str, nargs='+', metavar='TABLE:FASTA')

    args = p.parse_args()

    input_paths = {}
    for arg in args.hits:
        genome_id, hmmer_path, fasta_path = arg.strip().split(':')
        input_paths[genome_id] = (hmmer_path, fasta_path)

    log("Identifying single copy genes in each HMMER output table.")
    single_copy_sets = []
    hmmer_tables = {}
    for genome_id in input_paths:
        hmmer_path, _ = input_paths[genome_id]
        table = pd.read_table(hmmer_path)
        hmmer_tables[genome_id] = table
        single_copies = table.groupby(['model_id']).score.count()[lambda x: x==1].index
        single_copy_sets.append(set(single_copies))

    log('Identifying shared single-copy genes.')
    shared_single_copy_genes = single_copy_sets[0]
    for gene_set in single_copy_sets[1:]:
        shared_single_copy_genes &= gene_set

    # TODO: Check that different single copy genes don't have any duplicate hits?
    # Some ORFs do have duplicate hits, but those appear to be multiple,
    # non-overlapping domains in the same protein.
    all_genes = {}
    log('Outputting list of shared single-copy genes.')
    with open(args.list_output, 'w') as list_handle:
        for gene in shared_single_copy_genes:
            all_hits = {}
            for genome_id in hmmer_tables:
                hit = hmmer_tables[genome_id][lambda x: x.model_id == gene].orf_id
                assert hit.shape == (1,), 'Dimensionality of results are wrong ({})'.format(hit.shape)
                all_hits[genome_id] = hit.values[0]
            all_genes[gene] = all_hits
            print(gene, file=list_handle)

    log('Loading FASTA files.')
    fasta_indices = {}
    for genome_id in input_paths:
        _, fasta_path = input_paths[genome_id]
        fasta_indices[genome_id] = index_seqs(fasta_path, 'fasta')

    log('Writing gene outputs.')
    with open(args.tsv_output, 'w') as tsv_handle:
        for gene in all_genes:
            with open(args.fasta_output_template.format(gene), 'w') as fasta_handle:
                for genome_id in all_genes[gene]:
                    orf_id = all_genes[gene][genome_id]
                    rec = fasta_indices[genome_id][orf_id]
                    seq_id = rec.id
                    seq_seq = str(rec.seq)
                    print(f'>{seq_id}\n{seq_seq}', file=fasta_handle)
                    print(genome_id, gene, orf_id, sep='\t', file=tsv_handle)
