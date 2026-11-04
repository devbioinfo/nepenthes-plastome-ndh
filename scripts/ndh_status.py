#!/usr/bin/env python3
"""Score the 11 plastid ndh loci (ndhA-ndhK) in each plastome and assign retention classes.

Annotation-based scoring (Section 2.5 criteria, applied automatically):
  F  functional   : CDS present, not flagged /pseudo, length a multiple of 3,
                    starts with ATG/ACG/GTG, ends in a stop codon, no internal stop,
                    and >= 80% of the reference length (--min-frac)
  P  pseudogene   : annotated (CDS or gene feature) but fails any test above
  -  absent       : no feature for the locus

Classes: ndh-complete (11 F), ndh-partial (1-10 F), ndh-null (0 F).

Automatic scores are a starting point; the published matrix (Figure 3, Table 1)
was curated manually, and differences should be inspected by hand.

Usage:
    python scripts/ndh_status.py data/genbank/*.gb --out results/ndh_status.csv
"""
import argparse, csv, os, sys
sys.path.insert(0, os.path.dirname(__file__))
from genbank_utils import read_genbank, cds_records

NDH = ['ndhA', 'ndhB', 'ndhC', 'ndhD', 'ndhE', 'ndhF', 'ndhG', 'ndhH', 'ndhI', 'ndhJ', 'ndhK']
# approximate angiosperm reference CDS lengths (bp, incl. stop), used only for truncation checks
REF_LEN = {'ndhA': 1092, 'ndhB': 1533, 'ndhC': 363, 'ndhD': 1503, 'ndhE': 306, 'ndhF': 2241,
           'ndhG': 531, 'ndhH': 1182, 'ndhI': 546, 'ndhJ': 477, 'ndhK': 678}
STOPS = {'TAA', 'TAG', 'TGA'}


def orf_ok(nt, gene, min_frac):
    nt = nt.upper()
    if len(nt) % 3 or len(nt) < 3 * 30:
        return False
    if nt[:3] not in ('ATG', 'ACG', 'GTG'):        # ACG start is RNA-edited in many plastid genes
        return False
    cod = [nt[i:i + 3] for i in range(0, len(nt), 3)]
    if cod[-1] not in STOPS or any(c in STOPS for c in cod[:-1]):
        return False
    return len(nt) >= min_frac * REF_LEN[gene]


def score_record(rec, min_frac):
    status = {g: '-' for g in NDH}
    for f in rec['features']:
        if f['type'] in ('gene', 'CDS'):
            g = (f['qualifiers'].get('gene') or [''])[0]
            g = g[0:3].lower() + g[3:] if g.lower().startswith('ndh') else g
            if g in status and status[g] == '-':
                status[g] = 'P'
    for gene, nt, pseudo in cds_records(rec):
        g = gene[0:3].lower() + gene[3:] if gene.lower().startswith('ndh') else gene
        if g in status and not pseudo and orf_ok(nt, g, min_frac):
            status[g] = 'F'                            # any intact copy (e.g. either IR copy of ndhB) counts
    return status


def klass(n_func):
    return 'ndh-complete' if n_func == 11 else ('ndh-null' if n_func == 0 else 'ndh-partial')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('genbank', nargs='+')
    ap.add_argument('--out', default='results/ndh_status.csv')
    ap.add_argument('--min-frac', type=float, default=0.8)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or '.', exist_ok=True)
    rows = []
    for path in a.genbank:
        for rec in read_genbank(path):
            s = score_record(rec, a.min_frac)
            nF = sum(v == 'F' for v in s.values()); nP = sum(v == 'P' for v in s.values())
            nA = sum(v == '-' for v in s.values())
            assert nF + nP + nA == 11
            rows.append([rec['organism'], rec['accession']] + [s[g] for g in NDH] + [nF, nP, nA, klass(nF)])
    with open(a.out, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['species', 'accession'] + NDH + ['n_functional', 'n_pseudogene', 'n_absent', 'class'])
        w.writerows(sorted(rows))
    print('wrote', a.out, len(rows), 'plastomes')


if __name__ == '__main__':
    main()
