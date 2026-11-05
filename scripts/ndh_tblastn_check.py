#!/usr/bin/env python3
"""Search each plastome for ndh loci the annotation lacks (tblastn, N. khasiana proteins as queries).

For every locus scored '-' (not annotated), the genome is searched with tblastn.
HSPs with e-value <= 1e-10 and identity >= 50% (to exclude paralogous ndhB/D/F hits)
are merged; query coverage >= 0.2 -> pseudogene fragment ('P'); coverage >= 0.8 -> flagged
for manual inspection as a possibly unannotated intact gene; otherwise absent.

Usage: python scripts/ndh_tblastn_check.py results/ndh_status_auto.csv data/genbank --ref MH923233
"""
import csv, os, subprocess, sys, tempfile
sys.path.insert(0, os.path.dirname(__file__))
from genbank_utils import read_genbank, extract
from ndh_status import NDH, gname

CODE = {}
B = 'TCAG'; AA = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
for i, a in enumerate(B):
    for j, b in enumerate(B):
        for k, c in enumerate(B):
            CODE[a + b + c] = AA[16 * i + 4 * j + k]
tr = lambda s: ''.join(CODE.get(s[i:i + 3], 'X') for i in range(0, len(s) - 2, 3)).rstrip('*')


def main(table, gbdir, ref):
    refrec = next(read_genbank(os.path.join(gbdir, ref + '.gb')))
    prot = {}
    for f in refrec['features']:
        g = gname(f['qualifiers'])
        if f['type'] == 'CDS' and g in NDH and g not in prot:
            prot[g] = tr(extract(refrec['seq'], f['location']).upper())
    tmp = tempfile.mkdtemp()
    rows = list(csv.DictReader(open(table)))
    out = []
    for r in rows:
        acc = r['accession'].split('.')[0]
        rec = next(read_genbank(os.path.join(gbdir, acc + '.gb')))
        gfa = os.path.join(tmp, 'g.fa'); open(gfa, 'w').write('>g\n' + rec['seq'] + '\n')
        for g in NDH:
            if r[g] != '-':
                continue
            qfa = os.path.join(tmp, 'q.fa'); open(qfa, 'w').write(f'>{g}\n{prot[g]}\n')
            res = subprocess.run(['tblastn', '-query', qfa, '-subject', gfa, '-evalue', '1e-10',
                                  '-db_gencode', '11', '-outfmt', '6 qstart qend pident evalue sstart send'],
                                 capture_output=True, text=True).stdout
            cov = set()
            for line in res.strip().splitlines():
                qs, qe, pid = line.split()[:3]
                if float(pid) >= 50:
                    cov.update(range(int(qs), int(qe) + 1))
            c = len(cov) / len(prot[g])
            call = 'CHECK' if c >= 0.8 else ('P' if c >= 0.2 else '-')
            out.append([r['species'], g, f'{c:.2f}', call])
    w = csv.writer(sys.stdout); w.writerow(['species', 'gene', 'query_coverage', 'call'])
    w.writerows(out)


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument('table'); ap.add_argument('gbdir'); ap.add_argument('--ref', default='MH923233')
    a = ap.parse_args(); main(a.table, a.gbdir, a.ref)
