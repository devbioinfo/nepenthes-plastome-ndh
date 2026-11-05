#!/usr/bin/env python3
"""Effective number of codons (ENC; Wright 1990) and GC3s for each plastome (Figure 5c).

Uses the concatenated intact protein-coding sequences (not /pseudo, length a multiple of 3,
no internal stop). GC3s = GC content at synonymous third positions (Met, Trp and stops excluded).
ENC = 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6, where Fk is the mean homozygosity of k-fold amino acids;
if no 3-fold value exists, F3 is the mean of F2 and F4. Expected ENC under mutation pressure alone:
ENC = 2 + s + 29/(s^2 + (1 - s)^2), with s = GC3s.

Usage: python scripts/enc_gc3.py data/genbank/*.gb --out results/codon/enc_gc3.csv
"""
import argparse, csv, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from genbank_utils import read_genbank, cds_records
from rscu_pca import CODE

STOPS = {'TAA', 'TAG', 'TGA'}


def intact(nt):
    if len(nt) % 3: return False
    c = [nt[i:i + 3] for i in range(0, len(nt), 3)]
    return not any(x in STOPS for x in c[:-1])


def enc_gc3(seqs):
    counts = {c: 0 for c in CODE}
    for nt in seqs:
        for i in range(0, len(nt) - 2, 3):
            c = nt[i:i + 3]
            if c in counts: counts[c] += 1
    fam = {}
    for c, a in CODE.items():
        if a not in '*MW': fam.setdefault(a, []).append(c)
    F = {}
    for a, cods in fam.items():
        n = sum(counts[c] for c in cods)
        if n < 2: continue
        p = np.array([counts[c] / n for c in cods])
        F.setdefault(len(cods), []).append((n * (p ** 2).sum() - 1) / (n - 1))
    Fm = {k: np.mean(v) for k, v in F.items()}
    if 3 not in Fm: Fm[3] = (Fm[2] + Fm[4]) / 2
    enc = 2 + 9 / Fm[2] + 1 / Fm[3] + 5 / Fm[4] + 3 / Fm[6]
    syn = [c for c in CODE if CODE[c] not in '*MW']
    tot = sum(counts[c] for c in syn); gc = sum(counts[c] for c in syn if c[2] in 'GC')
    return min(enc, 61.0), gc / tot


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('genbank', nargs='+'); ap.add_argument('--out', default='results/codon/enc_gc3.csv')
    a = ap.parse_args()
    rows = []
    for p in a.genbank:
        for rec in read_genbank(p):
            seqs = [nt.upper() for g, nt, ps in cds_records(rec) if not ps and intact(nt.upper())]
            e, s = enc_gc3(seqs)
            exp = 2 + s + 29 / (s ** 2 + (1 - s) ** 2)
            rows.append([rec['organism'], f'{e:.2f}', f'{s:.4f}', f'{exp:.2f}', f'{(exp - e) / exp:.4f}'])
    with open(a.out, 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['species', 'ENC', 'GC3s', 'ENC_expected', 'ENC_ratio']); w.writerows(sorted(rows))
    print('wrote', a.out)


if __name__ == '__main__':
    main()
