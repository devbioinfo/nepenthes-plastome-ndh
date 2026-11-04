#!/usr/bin/env python3
"""Check a curated ndh status matrix: every row must account for all 11 loci.

Input CSV columns: species, ndhA ... ndhK with values F (functional), P (pseudogene), - (absent).
Prints per-row counts and the ndh retention class, and exits non-zero on any bad row.

Usage:
    python scripts/check_ndh_matrix.py data/ndh_status_curated.csv
"""
import csv, sys

NDH = ['ndhA', 'ndhB', 'ndhC', 'ndhD', 'ndhE', 'ndhF', 'ndhG', 'ndhH', 'ndhI', 'ndhJ', 'ndhK']


def main(path):
    bad = 0
    with open(path) as fh:
        for r in csv.DictReader(fh):
            vals = [r.get(g, '').strip() for g in NDH]
            if any(v not in ('F', 'P', '-') for v in vals):
                print('INCOMPLETE', r['species'], vals); bad += 1; continue
            nF, nP, nA = vals.count('F'), vals.count('P'), vals.count('-')
            cls = 'ndh-complete' if nF == 11 else ('ndh-null' if nF == 0 else 'ndh-partial')
            absent = ', '.join(g for g, v in zip(NDH, vals) if v == '-') or 'none'
            print(f"{r['species']:<40} F={nF:2d} P={nP:2d} absent={nA:2d} ({absent})  {cls}")
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main(sys.argv[1])
