#!/usr/bin/env python3
"""Print why each ndh locus was scored as it was (for manual curation)."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from genbank_utils import read_genbank, extract
from ndh_status import NDH, REF_LEN, STOPS

def diag(nt, g):
    nt = nt.upper(); cod = [nt[i:i+3] for i in range(0, len(nt) - len(nt) % 3, 3)]
    istops = [i for i, c in enumerate(cod[:-1]) if c in STOPS]
    return (f"len={len(nt)} ({len(nt)/REF_LEN[g]:.2f} ref) mod3={len(nt)%3} start={nt[:3]} end={cod[-1] if cod else ''} "
            f"intStops={len(istops)}{' first@codon '+str(istops[0]) if istops else ''}")

for path in sys.argv[1:]:
    for rec in read_genbank(path):
        print('==', rec['organism'], rec['accession'])
        for f in rec['features']:
            q = f['qualifiers']; g = (q.get('gene') or q.get('product') or [''])[0]
            if f['type'] in ('CDS', 'gene') and 'ndh' in str(q).lower():
                gg = g[:3].lower() + g[3:] if g.lower().startswith('ndh') else g
                extra = 'PSEUDO ' if ('pseudo' in q or 'pseudogene' in q) else ''
                d = diag(extract(rec['seq'], f['location']), gg) if (f['type']=='CDS' and gg in REF_LEN) else ''
                print(f"  {f['type']:4} {g:12} {extra}{f['location'][:60]:60} {d}")
