"""Minimal GenBank flat-file reader for plastome analyses (no Biopython needed).

Reads the sequence and CDS/gene features of one or more records, and extracts
CDS nucleotide sequences, handling complement(), join() and order() locations.
"""
import re

COMP = str.maketrans('ACGTRYKMSWBDHVNacgtrykmswbdhvn', 'TGCAYRMKSWVHDBNtgcayrmkswvhdbn')


def revcomp(s):
    return s.translate(COMP)[::-1]


def _segments(loc):
    """Parse a GenBank location into [(start0, end, strand)] in transcription order."""
    loc = loc.replace(' ', '')
    if loc.startswith('complement(') and loc.endswith(')'):
        return [(a, b, -s) for a, b, s in reversed(_segments(loc[11:-1]))]
    m = re.match(r'^(join|order)\((.*)\)$', loc)
    if m:
        out = []
        for piece in _split_top(m.group(2)):
            out.extend(_segments(piece))
        return out
    loc = loc.replace('<', '').replace('>', '')
    if '..' in loc:
        a, b = loc.split('..')
        return [(int(a) - 1, int(b), 1)]
    n = int(re.sub(r'\D', '', loc))
    return [(n - 1, n, 1)]


def _split_top(s):
    out, depth, cur = [], 0, ''
    for ch in s:
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth -= 1
        if ch == ',' and depth == 0:
            out.append(cur); cur = ''
        else:
            cur += ch
    if cur:
        out.append(cur)
    return out


def extract(seq, loc):
    """Spliced sequence of a feature (handles complement, join, trans-splicing)."""
    out = ''
    for a, b, s in _segments(loc):
        frag = seq[a:b]
        out += revcomp(frag) if s == -1 else frag
    return out


def read_genbank(path):
    """Yield dicts: {accession, organism, seq, features:[{type, location, qualifiers}]}."""
    with open(path) as fh:
        text = fh.read()
    for rec in text.split('\n//'):
        if 'LOCUS' not in rec:
            continue
        acc = re.search(r'^VERSION\s+(\S+)', rec, re.M) or re.search(r'^LOCUS\s+(\S+)', rec, re.M)
        org = re.search(r'^\s+ORGANISM\s+(.+)$', rec, re.M)
        seq = ''
        if 'ORIGIN' in rec:
            seq = re.sub(r'[^A-Za-z]', '', rec.split('ORIGIN', 1)[1]).upper()
        feats = []
        if 'FEATURES' in rec:
            block = rec.split('FEATURES', 1)[1].split('ORIGIN', 1)[0]
            cur = None
            for line in block.splitlines()[1:]:
                if len(line) > 21 and line[5] != ' ':
                    cur = {'type': line[5:21].strip(), 'location': line[21:].strip(), 'qualifiers': {}}
                    feats.append(cur); lastq = None
                elif cur is not None and line.startswith(' ' * 21):
                    body = line[21:]
                    if body.startswith('/'):
                        if '=' in body:
                            k, v = body[1:].split('=', 1)
                        else:
                            k, v = body[1:], True
                        cur['qualifiers'].setdefault(k, []).append(v.strip('"') if isinstance(v, str) else v)
                        lastq = k
                    elif lastq is None:
                        cur['location'] += body.strip()
                    else:
                        vals = cur['qualifiers'][lastq]
                        if isinstance(vals[-1], str):
                            vals[-1] = (vals[-1] + ' ' + body.strip()).strip('"')
        yield {'accession': acc.group(1) if acc else '', 'organism': org.group(1).strip() if org else '',
               'seq': seq, 'features': feats}


def cds_records(rec):
    """Yield (gene, nucleotide_seq, is_pseudo) for each CDS feature."""
    for f in rec['features']:
        if f['type'] != 'CDS':
            continue
        q = f['qualifiers']
        gene = (q.get('gene') or q.get('locus_tag') or ['?'])[0]
        pseudo = 'pseudo' in q or 'pseudogene' in q
        try:
            nt = extract(rec['seq'], f['location'])
        except Exception:
            continue
        yield gene, nt, pseudo
