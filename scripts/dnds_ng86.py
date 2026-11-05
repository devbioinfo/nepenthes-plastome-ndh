#!/usr/bin/env python3
"""Pairwise dN, dS and dN/dS (Nei & Gojobori 1986, Jukes-Cantor corrected) against N. khasiana.

Genes: 11 core photosynthetic genes (atpA, atpB, petA, petB, psaA, psaB, psbA, psbB, psbC,
psbD, rbcL) and the 11 plastid ndh genes. For each gene, intact copies (no internal stop,
length a multiple of 3, not /pseudo; for ndh genes only loci scored functional in the curated
matrix) are translated, aligned with MAFFT L-INS-i and back-translated. Each species is then
compared with N. khasiana after removing codons with a gap in either sequence.
Comparisons covering < 70% of the reference codons or with < 60% amino-acid identity are
discarded as partial or misannotated copies. Pairs with dS = 0 give an undefined ratio and
are excluded from ratio summaries; pairs with dS > 2 are flagged as saturated.

Usage: python scripts/dnds_ng86.py data/genbank/*.gb --curated results/ndh_status_curated.csv --outdir results/dnds
"""
import argparse, csv, os, subprocess, sys, warnings
sys.path.insert(0, os.path.dirname(__file__))
from genbank_utils import read_genbank, extract
from build_cds_alignment import CODE, STOPS, norm, taxon_name, intact
warnings.filterwarnings('ignore')
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds

CORE = ['atpA', 'atpB', 'petA', 'petB', 'psaA', 'psaB', 'psbA', 'psbB', 'psbC', 'psbD', 'rbcL']
NDH = ['ndhA', 'ndhB', 'ndhC', 'ndhD', 'ndhE', 'ndhF', 'ndhG', 'ndhH', 'ndhI', 'ndhJ', 'ndhK']
REF = 'Nepenthes_khasiana'


def codon_align(seqs, path, threads='2'):
    faa = path + '.faa'
    with open(faa, 'w') as fh:
        for t, nt in seqs.items():
            fh.write(f'>{t}\n' + ''.join(CODE.get(nt[i:i + 3], 'X') for i in range(0, len(nt), 3)) + '\n')
    aln = subprocess.run(['mafft', '--localpair', '--maxiterate', '1000', '--thread', threads, '--quiet', faa],
                         capture_output=True, text=True, check=True).stdout
    out = {}
    for block in aln.split('>')[1:]:
        name, *ls = block.splitlines(); p = ''.join(ls); nt = seqs[name.strip()]; k = 0; s = []
        for ch in p:
            if ch == '-': s.append('---')
            else: s.append(nt[k:k + 3]); k += 3
        out[name.strip()] = s
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('genbank', nargs='+')
    ap.add_argument('--curated', default='results/ndh_status_curated.csv')
    ap.add_argument('--outdir', default='results/dnds')
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    func = {}
    for r in csv.DictReader(open(a.curated)):
        func[taxon_name(r['species'])] = {g for g in NDH if r[g] == 'F'}
    genes = {}
    for p in a.genbank:
        for rec in read_genbank(p):
            t = taxon_name(rec['organism'])
            for f in rec['features']:
                q = f['qualifiers']
                if f['type'] != 'CDS' or 'pseudo' in q:
                    continue
                g = norm((q.get('gene') or q.get('product') or [''])[0])
                if g not in CORE + NDH or t in genes.get(g, {}):
                    continue
                if g in NDH and g not in func.get(t, set()):
                    continue
                nt = extract(rec['seq'], f['location']).upper()
                if intact(nt):
                    genes.setdefault(g, {})[t] = nt[:-3] if nt[-3:] in STOPS else nt
    rows = []
    for g in CORE + NDH:
        seqs = genes.get(g, {})
        if REF not in seqs or len(seqs) < 2:
            continue
        aln = codon_align(seqs, os.path.join(a.outdir, g))
        ref = aln[REF]
        for t, s in aln.items():
            if t == REF:
                continue
            keep = [i for i in range(len(ref)) if ref[i] != '---' and s[i] != '---'
                    and 'N' not in ref[i] + s[i] and ref[i] not in STOPS and s[i] not in STOPS]
            nref = sum(c != '---' for c in ref)
            if len(keep) < 0.7 * nref:
                continue                     # partial or misannotated copy
            aai = sum(CODE.get(ref[i], 'X') == CODE.get(s[i], 'Y') for i in keep) / len(keep)
            if aai < 0.6:
                continue                     # likely misannotation or misalignment
            c1 = CodonSeq(''.join(ref[i] for i in keep)); c2 = CodonSeq(''.join(s[i] for i in keep))
            try:
                dn, ds = cal_dn_ds(c1, c2, method='NG86')
            except Exception:
                continue
            ratio = dn / ds if ds > 0 else ''
            rows.append([g, 'ndh' if g in NDH else 'core', t, len(keep), f'{dn:.5f}', f'{ds:.5f}',
                         f'{ratio:.4f}' if ratio != '' else '', 'saturated' if ds > 2 else ''])
        print(g, len(seqs), 'sequences', flush=True)
    with open(os.path.join(a.outdir, 'dnds_vs_khasiana.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['gene', 'class', 'species', 'codons', 'dN', 'dS', 'dN_dS', 'flag'])
        w.writerows(rows)


if __name__ == '__main__':
    main()
