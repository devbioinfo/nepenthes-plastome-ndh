#!/usr/bin/env python3
"""Download the 50 plastome GenBank records listed in data/plastome_accessions.csv.

Usage:
    python scripts/fetch_plastomes.py --email you@example.org [--out data/genbank]

Uses NCBI E-utilities (efetch). Please supply your own e-mail as NCBI requests.
"""
import argparse, csv, os, time, urllib.parse, urllib.request

EFETCH = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--accessions', default='data/plastome_accessions.csv')
    ap.add_argument('--out', default='data/genbank')
    ap.add_argument('--email', required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    with open(a.accessions) as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        acc = r['genbank_accession'].split()[0]          # e.g. "MH923233 (NC_051455)" -> MH923233
        dest = os.path.join(a.out, acc + '.gb')
        if os.path.exists(dest) and os.path.getsize(dest) > 1000:
            continue
        q = urllib.parse.urlencode({'db': 'nuccore', 'id': acc, 'rettype': 'gbwithparts',
                                    'retmode': 'text', 'email': a.email, 'tool': 'carnivore_plastome'})
        with urllib.request.urlopen(EFETCH + '?' + q) as resp:
            data = resp.read().decode()
        with open(dest, 'w') as fh:
            fh.write(data)
        print('fetched', acc, r['species'])
        time.sleep(0.4)                                   # stay under 3 requests/s


if __name__ == '__main__':
    main()
