#!/usr/bin/env python3
"""RSCU from concatenated plastid CDS and PCA (independent cross-check of Figure 5).

RSCU(codon) = observed count / (total count for its amino acid / number of synonymous codons).
Stop codons and the single-codon amino acids Met (ATG) and Trp (TGG) are excluded by default
(59 codons); --keep-trp keeps TGG to reproduce a 60-codon matrix.
PCA is computed on the centred and scaled species x codon matrix (numpy SVD).

Usage:
    python scripts/rscu_pca.py data/genbank/*.gb --exclude-outgroup Amborella --outdir results/codon
"""
import argparse, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from genbank_utils import read_genbank, cds_records

BASES = 'TCAG'
AA = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODE = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate(BASES) for j, b in enumerate(BASES)
        for k, c in enumerate(BASES)}


def rscu(seqs, keep_trp=False):
    counts = {c: 0 for c in CODE}
    for nt in seqs:
        nt = nt.upper()
        for i in range(0, len(nt) - len(nt) % 3, 3):
            c = nt[i:i + 3]
            if c in counts:
                counts[c] += 1
    syn = {}
    for c, a in CODE.items():
        syn.setdefault(a, []).append(c)
    out = {}
    for a, cods in syn.items():
        if a == '*' or a == 'M' or (a == 'W' and not keep_trp):
            continue
        tot = sum(counts[c] for c in cods)
        for c in cods:
            out[c] = counts[c] / (tot / len(cods)) if tot else np.nan
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('genbank', nargs='+')
    ap.add_argument('--outdir', default='results/codon')
    ap.add_argument('--keep-trp', action='store_true')
    ap.add_argument('--exclude-outgroup', default='Amborella')
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    table = {}
    for path in a.genbank:
        for rec in read_genbank(path):
            if a.exclude_outgroup and rec['organism'].startswith(a.exclude_outgroup):
                continue
            seqs = [nt for g, nt, pseudo in cds_records(rec) if not pseudo and len(nt) % 3 == 0]
            table[rec['organism']] = rscu(seqs, a.keep_trp)
    df = pd.DataFrame(table).T.sort_index()
    df = df.dropna(axis=0, how='all')            # plastomes with no usable CDS
    df.to_csv(os.path.join(a.outdir, 'rscu.csv'))
    X = df.dropna(axis=1).values
    sd = X.std(axis=0, ddof=1); keep = sd > 0
    Z = (X[:, keep] - X[:, keep].mean(axis=0)) / sd[keep]
    U, S, Vt = np.linalg.svd(Z, full_matrices=False)
    var = S ** 2 / (S ** 2).sum()
    scores = pd.DataFrame(U * S, index=df.index, columns=[f'PC{i+1}' for i in range(len(S))])
    scores.iloc[:, :5].to_csv(os.path.join(a.outdir, 'pca_scores.csv'))
    pd.Series(var, index=scores.columns, name='variance_explained').to_csv(os.path.join(a.outdir, 'pca_variance.csv'))
    print(f'{df.shape[0]} plastomes x {df.shape[1]} codons; PC1 {var[0]:.1%}, PC2 {var[1]:.1%}')


if __name__ == '__main__':
    main()
