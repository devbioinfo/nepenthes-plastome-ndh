#!/usr/bin/env python3
"""Combine annotation-based scores with sequence evidence into the curated ndh matrix.

Rules (Methods 2.5):
 1. annotation-based score (ndh_status.py);
 2. loci with no annotation: tblastn fragment (coverage 0.2-0.8) -> P; >= 0.8 resolved by ORF test;
 3. loci annotated as gene only (no CDS, not /pseudo): resolved by the same ORF test;
 4. loci flagged /pseudo by the submitters stay P.
"""
import csv, sys
NDH = ['ndhA','ndhB','ndhC','ndhD','ndhE','ndhF','ndhG','ndhH','ndhI','ndhJ','ndhK']
auto = list(csv.DictReader(open('results/ndh_status_auto.csv')))
reasons = {(r['species'], r['gene']): r['reason'] for r in csv.DictReader(open('results/ndh_status_auto_reasons.csv'))}
tb = {(r['species'], r['gene']): r['call'] for r in csv.DictReader(open('results/ndh_tblastn_absent_loci.csv'))}
res = {(r['species'], r['gene']): r['call'] for r in csv.DictReader(open('results/ndh_unannotated_resolved.csv'))}
pre = {(r['species'], r['gene']): r['call'] for r in csv.DictReader(open('results/ndh_pseudogene_recheck.csv'))}
changes = []
out = []
for r in auto:
    sp = r['species']; row = {g: r[g] for g in NDH}
    for g in NDH:
        k = (sp, g); old = row[g]
        if old == '-' and tb.get(k) == 'P':
            row[g] = 'P'
        elif old == '-' and tb.get(k) == 'CHECK':
            row[g] = res[k]
        elif old == 'P' and pre.get(k) == 'F' and reasons[k] != 'flagged /pseudo':
            # gene feature only; check it was not /pseudo in the record
            row[g] = 'F'
        if row[g] != old:
            changes.append([sp, g, old, row[g]])
    nF = sum(v == 'F' for v in row.values()); nP = sum(v == 'P' for v in row.values()); nA = 11 - nF - nP
    cls = 'ndh-complete' if nF == 11 else ('ndh-null' if nF == 0 else 'ndh-partial')
    out.append([sp, r['accession']] + [row[g] for g in NDH] + [nF, nP, nA, cls])
with open('results/ndh_status_curated.csv', 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['species', 'accession'] + NDH + ['n_functional', 'n_pseudogene', 'n_absent', 'class']); w.writerows(out)
with open('results/ndh_curation_changes.csv', 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['species', 'gene', 'annotation_score', 'curated']); w.writerows(changes)
print(len(changes), 'changes')
