#!/usr/bin/env python3
"""Alignment-free 15-mer distances between plastomes and a neighbour-joining tree (Figure 3).

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


def nj(D, names):
    """Saitou & Nei neighbour-joining; returns an unrooted Newick string."""
    D = D.copy().astype(float)
    nodes = [n.replace(' ', '_') for n in names]
    while len(nodes) > 3:
        n = len(nodes)
        r = D.sum(axis=1)
        Q = (n - 2) * D - r[:, None] - r[None, :]
        np.fill_diagonal(Q, np.inf)
        i, j = np.unravel_index(np.argmin(Q), Q.shape)
        if i > j: i, j = j, i
        li = 0.5 * D[i, j] + (r[i] - r[j]) / (2 * (n - 2))
        lj = D[i, j] - li
        new = f'({nodes[i]}:{max(li, 0):.6f},{nodes[j]}:{max(lj, 0):.6f})'
        du = 0.5 * (D[i] + D[j] - D[i, j])
        keep = [k for k in range(n) if k not in (i, j)]
        D2 = np.zeros((len(keep) + 1, len(keep) + 1))
        D2[:-1, :-1] = D[np.ix_(keep, keep)]
        D2[-1, :-1] = D2[:-1, -1] = du[keep]
        D = D2
        nodes = [nodes[k] for k in keep] + [new]
    a, b, c = nodes
    la = 0.5 * (D[0, 1] + D[0, 2] - D[1, 2]); lb = D[0, 1] - la; lc = D[0, 2] - la
    return f'({a}:{max(la,0):.6f},{b}:{max(lb,0):.6f},{c}:{max(lc,0):.6f});'


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
    with open(os.path.join(a.outdir, 'kmer_nj.nwk'), 'w') as fh:
        fh.write(nj(D, names) + '\n')
    print('wrote distances, NJ and UPGMA trees for', n, 'plastomes')


if __name__ == '__main__':
    main()
