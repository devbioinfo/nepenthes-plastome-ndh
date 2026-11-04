#!/usr/bin/env python3
"""Alignment-free 15-mer distances between plastomes and a UPGMA tree (cross-check of Figure 3).

Distance = Mash-style estimate from the Jaccard index J of canonical k-mer sets:
    D = -1/k * ln(2J / (1 + J))

Usage:
    python scripts/kmer_distance.py data/genbank/*.gb --k 15 --outdir results/kmer
"""
import argparse, math, os, sys
import numpy as np
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform
sys.path.insert(0, os.path.dirname(__file__))
from genbank_utils import read_genbank, revcomp


def kmers(seq, k):
    s = set()
    seq = seq.upper()
    rc = revcomp(seq)
    n = len(seq)
    for i in range(n - k + 1):
        f = seq[i:i + k]
        if 'N' in f:
            continue
        r = rc[n - i - k:n - i]
        s.add(min(f, r))
    return s


def newick(node, names):
    # UPGMA node height = merge distance / 2; branch length = parent height - child height
    if node.is_leaf():
        return names[node.id].replace(' ', '_')
    h = node.dist / 2
    parts = []
    for ch in (node.get_left(), node.get_right()):
        parts.append(f'{newick(ch, names)}:{h - ch.dist / 2:.6f}')
    return '(' + ','.join(parts) + ')'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('genbank', nargs='+')
    ap.add_argument('--k', type=int, default=15)
    ap.add_argument('--outdir', default='results/kmer')
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    names, sets = [], []
    for path in a.genbank:
        for rec in read_genbank(path):
            names.append(rec['organism']); sets.append(kmers(rec['seq'], a.k))
    n = len(names)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            inter = len(sets[i] & sets[j]); union = len(sets[i] | sets[j])
            J = inter / union if union else 0
            d = 1.0 if J == 0 else -1.0 / a.k * math.log(2 * J / (1 + J))
            D[i, j] = D[j, i] = d
    with open(os.path.join(a.outdir, 'kmer_distances.tsv'), 'w') as fh:
        fh.write('\t' + '\t'.join(names) + '\n')
        for i in range(n):
            fh.write(names[i] + '\t' + '\t'.join(f'{x:.6f}' for x in D[i]) + '\n')
    Z = linkage(squareform(D, checks=False), method='average')
    with open(os.path.join(a.outdir, 'kmer_upgma.nwk'), 'w') as fh:
        fh.write(newick(to_tree(Z), names) + ';\n')
    print('wrote distances and UPGMA tree for', n, 'plastomes')


if __name__ == '__main__':
    main()
