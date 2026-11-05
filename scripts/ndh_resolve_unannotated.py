#!/usr/bin/env python3
"""Resolve unannotated ndh loci with near-full tblastn coverage (call = CHECK).

Intact ('F') if: tblastn HSPs (identity >= 50%) cover >= 80% of the N. khasiana protein,
no HSP contains a stop codon, and for single-exon genes an open reading frame
(ATG/ACG/GTG ... stop) spanning the hit is >= 80% of the reference length.
Otherwise pseudogene ('P'). ndhA and ndhB contain an intron, so for them the ORF test
is applied to the HSPs only (no stops, consistent frame per exon).
"""
import csv, os, subprocess, sys, tempfile
sys.path.insert(0, os.path.dirname(__file__))
from genbank_utils import read_genbank, extract, revcomp
from ndh_status import NDH, REF_LEN, STOPS, gname
from ndh_tblastn_check import tr

def orf_around(seq, s, e, strand):
    """ORF containing [s,e) (0-based, on given strand); returns length in nt or 0."""
    if strand == -1:
        seq = revcomp(seq); s, e = len(seq) - e, len(seq) - s
    frame = s % 3
    # walk upstream from s to a stop; last start codon seen is ORF start
    i = s; start = None
    while i - 3 >= 0:
        i -= 3; c = seq[i:i+3]
        if c in STOPS: break
        if c in ('ATG', 'ACG', 'GTG'): start = i
    if start is None:
        start = s if seq[s:s+3] in ('ATG', 'ACG', 'GTG') else None
    j = s
    while j + 3 <= len(seq):
        c = seq[j:j+3]
        if c in STOPS: break
        j += 3
    if start is None or j < e - 3:
        return 0
    return j + 3 - start

def main(check_csv, status_csv, gbdir, ref='MH923233'):
    refrec = next(read_genbank(os.path.join(gbdir, ref + '.gb')))
    prot = {}
    for f in refrec['features']:
        g = gname(f['qualifiers'])
        if f['type'] == 'CDS' and g in NDH and g not in prot:
            prot[g] = tr(extract(refrec['seq'], f['location']).upper())
    accs = {r['species']: r['accession'].split('.')[0] for r in csv.DictReader(open(status_csv))}
    tmp = tempfile.mkdtemp(); out = []
    for r in csv.DictReader(open(check_csv)):
        if r['call'] != 'CHECK':
            continue
        sp, g = r['species'], r['gene']
        rec = next(read_genbank(os.path.join(gbdir, accs[sp] + '.gb'))); seq = rec['seq'].upper()
        gfa = os.path.join(tmp, 'g.fa'); open(gfa, 'w').write('>g\n' + seq + '\n')
        qfa = os.path.join(tmp, 'q.fa'); open(qfa, 'w').write(f'>{g}\n{prot[g]}\n')
        res = subprocess.run(['tblastn', '-query', qfa, '-subject', gfa, '-evalue', '1e-10', '-db_gencode', '11',
                              '-outfmt', '6 qstart qend pident sstart send sframe sseq'], capture_output=True, text=True).stdout
        hsps = [l.split('\t') for l in res.strip().splitlines() if float(l.split('\t')[2]) >= 50]
        # keep the best locus only (IR copies give duplicate hits): group by subject region
        hsps.sort(key=lambda h: -(int(h[1]) - int(h[0])))
        if not hsps:
            out.append([sp, g, '0.00', 'P', 'no HSP with identity >= 50%']); continue
        top = hsps[0]; region = (min(int(top[3]), int(top[4])))
        loc_hsps = [h for h in hsps if abs(min(int(h[3]), int(h[4])) - region) < 4000]
        stops = sum(h[6].count('*') for h in loc_hsps)
        frames = {h[5] for h in loc_hsps}
        cov = set()
        for h in loc_hsps: cov.update(range(int(h[0]), int(h[1]) + 1))
        c = len(cov) / len(prot[g])
        orf = ''
        if g in ('ndhA', 'ndhB'):
            call = 'F' if (stops == 0 and c >= 0.8) else 'P'
            note = f'intron gene; HSP stops={stops}, frames={sorted(frames)}'
        else:
            s = min(min(int(h[3]), int(h[4])) for h in loc_hsps) - 1
            e = max(max(int(h[3]), int(h[4])) for h in loc_hsps)
            strand = 1 if int(top[5]) > 0 else -1
            L = orf_around(seq, s, e, strand) if len(frames) == 1 else 0
            call = 'F' if (stops == 0 and c >= 0.8 and L >= 0.8 * REF_LEN[g]) else 'P'
            note = f'HSP stops={stops}, frames={sorted(frames)}, ORF={L} nt ({L/REF_LEN[g]:.2f} ref)'
        out.append([sp, g, f'{c:.2f}', call, note])
    w = csv.writer(sys.stdout); w.writerow(['species', 'gene', 'coverage', 'call', 'note']); w.writerows(out)

if __name__ == '__main__':
    main(*sys.argv[1:4])
