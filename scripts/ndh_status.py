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
from genbank_utils import read_genbank, cds_records, extract, revcomp, _segments

NDH = ['ndhA', 'ndhB', 'ndhC', 'ndhD', 'ndhE', 'ndhF', 'ndhG', 'ndhH', 'ndhI', 'ndhJ', 'ndhK']
# approximate angiosperm reference CDS lengths (bp, incl. stop), used only for truncation checks
REF_LEN = {'ndhA': 1092, 'ndhB': 1533, 'ndhC': 363, 'ndhD': 1503, 'ndhE': 306, 'ndhF': 2241,
           'ndhG': 531, 'ndhH': 1182, 'ndhI': 546, 'ndhJ': 477, 'ndhK': 678}
STOPS = {'TAA', 'TAG', 'TGA'}


def gname(q):
    g = (q.get('gene') or q.get('product') or [''])[0].strip()
    return g[:3].lower() + g[3:] if g.lower().startswith('ndh') else g


def extend_start(seq, loc, nt):
    """Extend a CDS upstream, in frame, to the most upstream ATG/ACG/GTG before a stop.

    Plastid ndhD (and some other genes) start with ACG, edited to AUG; GenBank
    annotations sometimes begin at a downstream ATG instead.
    """
    segs = _segments(loc)
    a, b, strand = segs[0]                      # first exon in transcription order
    best = nt
    for k in range(1, 400):                     # up to 400 codons upstream
        if strand == 1:
            s0 = a - 3 * k
            if s0 < 0: break
            cod = seq[s0:s0 + 3]
        else:
            e0 = b + 3 * k
            if e0 > len(seq): break
            cod = revcomp(seq[e0 - 3:e0])
        cod = cod.upper()
        if cod in STOPS:
            break
        if cod in ('ATG', 'ACG', 'GTG'):
            if strand == 1:
                best = seq[s0:a] + nt
            else:
                best = revcomp(seq[b:e0]) + nt
    return best


def orf_ok(nt, gene, min_frac):
    nt = nt.upper()
    if len(nt) % 3 or len(nt) < 3 * 30:
        return False, 'frame'
    cod = [nt[i:i + 3] for i in range(0, len(nt), 3)]
    if cod[-1] not in STOPS:
        return False, 'no terminal stop'
    if any(c in STOPS for c in cod[:-1]):
        return False, 'internal stop'
    if len(nt) < min_frac * REF_LEN[gene]:
        return False, f'truncated ({len(nt)/REF_LEN[gene]:.2f} of reference)'
    return True, f'intact ({len(nt)/REF_LEN[gene]:.2f})'


def score_record(rec, min_frac, notes=None):
    status = {g: '-' for g in NDH}
    why = {g: 'not annotated' for g in NDH}
    for f in rec['features']:
        if f['type'] in ('gene', 'CDS'):
            g = gname(f['qualifiers'])
            if g in status:
                q = f['qualifiers']
                if 'pseudo' in q or 'pseudogene' in q:
                    if status[g] != 'F':
                        status[g] = 'P'; why[g] = 'flagged /pseudo'
                elif status[g] == '-':
                    status[g] = 'P'
                    why[g] = 'annotated pseudogene/gene only' if f['type'] == 'gene' else why[g]
    for f in rec['features']:
        if f['type'] != 'CDS':
            continue
        g = gname(f['qualifiers'])
        if g not in status or status[g] == 'F':
            continue
        q = f['qualifiers']
        if 'pseudo' in q or 'pseudogene' in q:
            why[g] = 'flagged /pseudo'; continue
        nt = extract(rec['seq'], f['location'])
        ok, reason = orf_ok(nt, g, min_frac)
        if not ok and 'truncated' in reason:
            nt2 = extend_start(rec['seq'], f['location'], nt)
            if len(nt2) > len(nt):
                ok2, reason2 = orf_ok(nt2, g, min_frac)
                if ok2:
                    ok, reason = True, reason2 + ', start extended upstream'
        if ok:
            status[g] = 'F'                            # any intact copy (e.g. either IR copy of ndhB) counts
        why[g] = reason
    if notes is not None:
        notes.update(why)
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
    reasons = []
    for path in a.genbank:
        for rec in read_genbank(path):
            notes = {}
            s = score_record(rec, a.min_frac, notes)
            for g in NDH:
                reasons.append([rec['organism'], g, s[g], notes[g]])
            nF = sum(v == 'F' for v in s.values()); nP = sum(v == 'P' for v in s.values())
            nA = sum(v == '-' for v in s.values())
            assert nF + nP + nA == 11
            rows.append([rec['organism'], rec['accession']] + [s[g] for g in NDH] + [nF, nP, nA, klass(nF)])
    with open(a.out, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['species', 'accession'] + NDH + ['n_functional', 'n_pseudogene', 'n_absent', 'class'])
        w.writerows(sorted(rows))
    with open(a.out.replace('.csv', '_reasons.csv'), 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['species', 'gene', 'status', 'reason']); w.writerows(reasons)
    print('wrote', a.out, len(rows), 'plastomes')


if __name__ == '__main__':
    main()
