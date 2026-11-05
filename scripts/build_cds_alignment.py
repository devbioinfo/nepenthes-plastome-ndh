#!/usr/bin/env python3
"""Build a concatenated codon alignment of plastid protein-coding genes for phylogenomics.

Steps
 1. From each GenBank record take one intact copy of every annotated CDS
    (not /pseudo, length a multiple of 3, no internal stop codon).
 2. Keep genes present in at least --min-taxa plastomes.
 3. Translate, align proteins with MAFFT (L-INS-i, --localpair --maxiterate 1000),
    and back-translate to a codon alignment.
 4. Drop codon columns with more than --max-gap gap fraction.
 5. Concatenate and write a NEXUS partition file (one partition per gene).

Usage:
  python scripts/build_cds_alignment.py data/genbank/*.gb --outdir results/phylogeny
"""
import argparse, os, subprocess, sys
sys.path.insert(0, os.path.dirname(__file__))
from genbank_utils import read_genbank, extract

B = 'TCAG'; AA = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODE = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate(B) for j, b in enumerate(B) for k, c in enumerate(B)}
STOPS = {'TAA', 'TAG', 'TGA'}


def norm(g):
    g = g.strip()
    return g[:3].lower() + g[3:] if len(g) > 3 else g.lower()


def taxon_name(org):
    return org.replace(' x ', '_x_').replace(' ', '_').replace('.', '')


def intact(nt):
    if len(nt) % 3 or len(nt) < 90:
        return False
    cod = [nt[i:i + 3] for i in range(0, len(nt), 3)]
    return not any(c in STOPS for c in cod[:-1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('genbank', nargs='+')
    ap.add_argument('--outdir', default='results/phylogeny')
    ap.add_argument('--min-taxa', type=int, default=40)
    ap.add_argument('--max-gap', type=float, default=0.5)
    ap.add_argument('--threads', default='2')
    a = ap.parse_args()
    os.makedirs(os.path.join(a.outdir, 'genes'), exist_ok=True)
    genes = {}
    taxa = []
    for path in a.genbank:
        for rec in read_genbank(path):
            t = taxon_name(rec['organism']); taxa.append(t)
            seen = set()
            for f in rec['features']:
                q = f['qualifiers']
                if f['type'] != 'CDS' or 'pseudo' in q or 'pseudogene' in q:
                    continue
                g = norm((q.get('gene') or q.get('product') or [''])[0])
                if not g or g in seen:
                    continue
                try:
                    nt = extract(rec['seq'], f['location']).upper()
                except Exception:
                    continue
                if intact(nt):
                    genes.setdefault(g, {})[t] = nt[:-3] if nt[-3:] in STOPS else nt
                    seen.add(g)
    kept = sorted(g for g, d in genes.items() if len(d) >= a.min_taxa)
    print(f'{len(genes)} genes found; {len(kept)} present in >= {a.min_taxa} of {len(taxa)} plastomes')
    concat = {t: '' for t in taxa}
    parts = []
    pos = 1
    log = []
    for g in kept:
        seqs = genes[g]
        faa = os.path.join(a.outdir, 'genes', g + '.faa')
        with open(faa, 'w') as fh:
            for t, nt in seqs.items():
                fh.write(f'>{t}\n' + ''.join(CODE.get(nt[i:i + 3], 'X') for i in range(0, len(nt), 3)) + '\n')
        aln = subprocess.run(['mafft', '--localpair', '--maxiterate', '1000', '--thread', a.threads, '--quiet', faa],
                             capture_output=True, text=True, check=True).stdout
        prot = {}
        for block in aln.split('>')[1:]:
            name, *ls = block.splitlines(); prot[name.strip()] = ''.join(ls)
        L = len(next(iter(prot.values())))
        codon = {}
        for t, p in prot.items():
            nt = seqs[t]; k = 0; s = []
            for ch in p:
                if ch == '-':
                    s.append('---')
                else:
                    s.append(nt[k:k + 3]); k += 3
            codon[t] = s
        keep = [i for i in range(L) if sum(codon[t][i] == '---' for t in codon) / len(taxa) <= a.max_gap
                and sum(codon[t][i] == '---' for t in codon) + (len(taxa) - len(codon)) <= a.max_gap * len(taxa)]
        for t in taxa:
            concat[t] += ''.join(codon[t][i] for i in keep) if t in codon else '---' * len(keep)
        n = 3 * len(keep)
        with open(os.path.join(a.outdir, 'genes', g + '.codon.fasta'), 'w') as fh:
            for t in codon:
                fh.write(f'>{t}\n' + ''.join(codon[t][i] for i in keep) + '\n')
        parts.append((g, pos, pos + n - 1)); pos += n
        log.append([g, len(seqs), L, len(keep)])
    with open(os.path.join(a.outdir, 'concatenated_cds.fasta'), 'w') as fh:
        for t in taxa:
            fh.write(f'>{t}\n{concat[t]}\n')
    with open(os.path.join(a.outdir, 'partitions.nex'), 'w') as fh:
        fh.write('#nexus\nbegin sets;\n')
        for g, s, e in parts:
            fh.write(f'  charset {g} = {s}-{e};\n')
        fh.write('end;\n')
    with open(os.path.join(a.outdir, 'genes_used.tsv'), 'w') as fh:
        fh.write('gene\tn_taxa\taligned_codons\tcodons_kept\n')
        for r in log:
            fh.write('\t'.join(map(str, r)) + '\n')
    print(f'concatenated alignment: {len(taxa)} taxa x {pos - 1} nt ({len(kept)} genes)')


if __name__ == '__main__':
    main()
