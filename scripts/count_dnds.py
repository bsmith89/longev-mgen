#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio.SeqIO import index, parse
from itertools import product
import sys
import numpy as np
import argparse

NUCLEOTIDES = {'A', 'C', 'G', 'T'}
TRANSLATION = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'$', 'TAG':'$',
    'TGC':'C', 'TGT':'C', 'TGA':'$', 'TGG':'W',
    }


def codons(seq):
    assert len(seq) % 3 == 0
    for i in range(0, len(seq), 3):
        yield seq[i:i+3]


def is_indel(codon):
    if '-' in codon:
        assert set(codon) == set('-'), "Out-of-frame indel found."
        return True
    else:
        return False

def count_codons(seqA, seqB):
    assert len(seqA) == len(seqB), "{} != {}".format(len(seqA), len(seqB))
    count_nsyn, count_syn, pos_nsyn, pos_syn, indels = 0, 0, 0, 0, 0
    for codonA, codonB in zip(codons(seqA), codons(seqB)):
        indelA = is_indel(codonA)
        indelB = is_indel(codonB)
        if indelA and indelB:
            continue
        elif indelA != indelB:
            # XOR each codon is an indel.
            indels += 1
            continue
        cn, cs = count_changes(codonA, codonB)
        pnA, psA = count_positions(codonA)
        pnB, psB = count_positions(codonB)
        count_syn += cs
        count_nsyn += cn
        pos_nsyn += (pnA + pnB)
        pos_syn += (psA + psB)
    return count_nsyn, count_syn, pos_nsyn, pos_syn, indels


def count_changes(codonA, codonB):
    """Count (non-)synonymous differences in a codon.
    """
    count_nsyn, count_syn = 0, 0
    if codonA == codonB:
        return 0, 0
    for i, (nA, nB) in enumerate(zip(codonA, codonB)):
        if nA == nB:
            continue
        codonB_prime = list(codonB)
        codonB_prime[i] = nA
        if is_same_aa(''.join(codonB_prime), codonB):
            count_syn += 1
        else:
            count_nsyn += 1
    return count_nsyn, count_syn

def is_same_aa(codonA, codonB):
    aaA = TRANSLATION[codonA]
    aaB = TRANSLATION[codonB]
    # print(codonA, aaA, codonB, aaB, file=sys.stderr)
    return aaA == aaB


def _count_positions(codon):
    """Count possible (non-)synonymous changes in a codon resulting from one mutation.
    """
    count_nsyn, count_syn = 0, 0
    for i, n in enumerate(codon):
        codon_prime = list(codon)
        for a in NUCLEOTIDES - {n}:
            codon_prime = list(codon)
            codon_prime[i] = a
            if is_same_aa(''.join(codon_prime), codon):
                count_syn += 1
            else:
                count_nsyn += 1
    return count_nsyn, count_syn


# Cache positions
POSITIONS = {}
for codon in TRANSLATION:
    POSITIONS[codon] = _count_positions(codon)


def count_positions(codon):
    return POSITIONS[codon]


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('align1', metavar='FASTA1')
    p.add_argument('align2', metavar='FASTA2', nargs='?')
    args = p.parse_args()

    if args.align2:
        rec_index1 = index(args.align1, 'fasta')
        rec_index2 = index(args.align2, 'fasta')
        comparisons = []
        for idA, idB in zip(rec_index1.keys(), rec_index2.keys()):
            comparisons.append((idA, rec_index1[idA],
                                idB, rec_index2[idB]))
    else:
        rec_index = index(args.align1, 'fasta')
        comparisons = []
        ids = list(rec_index.keys())
        for i, idA in enumerate(ids):
            for j, idB in enumerate(ids[i+1:]):
                comparisons.append((idA, rec_index[idA],
                                    idB, rec_index[idB]))

    for idA, recA, idB, recB in comparisons:
        seqA = recA.seq
        seqB = recB.seq
        cN, cS, pN, pS, ind = count_codons(seqA, seqB)
        if pN == 0:
            dN = float('nan')
        else:
            dN = cN/pN
        if pS == 0:
            dS = float('nan')
        else:
            dS = cS/pS
        if (dS == 0) or np.isnan(dN):
            dnds = float('nan')
        else:
            dnds = (cN/pN)/(cS/pS)
        print(recA.id, recB.id, cN, cS, pN, pS, dN, dS, ind, dnds, sep='\t')
