#!/usr/bin/env python3
"""Pairwise whole-plastome identity (ANIb-style, Goris et al. 2007) and hierarchical clustering (Figure 6).

Each plastome is cut into consecutive 1,020-nt fragments, which are searched with
blastn (task dc-megablast, e-value 1e-5) against every other plastome.
For each query-subject pair, the best hit of each fragment is kept if it covers >= 70%
of the fragment at >= 30% identity. ANI(query->subject) is the mean identity of those
hits; the reported value for a pair is the mean of both directions. The fraction of
fragments counted (alignment fraction) is also reported.

Usage: python scripts/plastome_ani.py data/genbank/*.gb --outdir results/ani --threads 2
"""
import argparse, os, subprocess, sys, tempfile
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from genbank_utils import read_genbank
from build_cds_alignment import taxon_name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('genbank', nargs='+')
    ap.add_argument('--outdir', default='results/ani')
    ap.add_argument('--threads', default='2')
    ap.add_argument('--frag', type=int, default=1020)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    tmp = tempfile.mkdtemp()
    recs = [r for p in a.genbank for r in read_genbank(p)]
    names = [taxon_name(r['organism']) for r in recs]
    db = os.path.join(tmp, 'all.fa')
    with open(db, 'w') as fh:
        for n, r in zip(names, recs):
            fh.write(f'>{n}\n{r["seq"]}\n')
    subprocess.run(['makeblastdb', '-in', db, '-dbtype', 'nucl', '-out', db], check=True, capture_output=True)
    n = len(names); idx = {x: i for i, x in enumerate(names)}
    ani = np.full((n, n), np.nan); af = np.zeros((n, n))
    for i, (q, r) in enumerate(zip(names, recs)):
        frags = os.path.join(tmp, 'q.fa'); nf = 0
        with open(frags, 'w') as fh:
            s = r['seq']
            for k in range(0, len(s) - a.frag + 1, a.frag):
                fh.write(f'>f{k}\n{s[k:k + a.frag]}\n'); nf += 1
        out = subprocess.run(['blastn', '-task', 'dc-megablast', '-query', frags, '-db', db, '-evalue', '1e-5',
                              '-max_target_seqs', str(n), '-max_hsps', '1', '-num_threads', a.threads,
                              '-outfmt', '6 qseqid sseqid pident length qstart qend'],
                             capture_output=True, text=True, check=True).stdout
        best = {}
        for line in out.splitlines():
            f, sub, pid, L, qs, qe = line.split('\t')
            cov = (abs(int(qe) - int(qs)) + 1) / a.frag
            if sub == q or cov < 0.7 or float(pid) < 30:
                continue
            key = (f, sub)
            if key not in best or float(pid) > best[key]:
                best[key] = float(pid)
        per = {}
        for (f, sub), pid in best.items():
            per.setdefault(sub, []).append(pid)
        for sub, vals in per.items():
            j = idx[sub]; ani[i, j] = np.mean(vals); af[i, j] = len(vals) / nf
        print(f'{i + 1}/{n} {q}', flush=True)
    sym = np.nanmean(np.stack([ani, ani.T]), axis=0)
    np.fill_diagonal(sym, 100.0)
    afs = (af + af.T) / 2; np.fill_diagonal(afs, 1.0)
    for fn, M in (('ani_matrix.tsv', sym), ('alignment_fraction.tsv', afs)):
        with open(os.path.join(a.outdir, fn), 'w') as fh:
            fh.write('\t' + '\t'.join(names) + '\n')
            for k in range(n):
                fh.write(names[k] + '\t' + '\t'.join('' if np.isnan(x) else f'{x:.3f}' for x in M[k]) + '\n')


if __name__ == '__main__':
    main()
