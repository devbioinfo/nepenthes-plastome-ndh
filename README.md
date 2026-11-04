# nepenthes-plastome-ndh

Data and code for:

> Konhar R, Dash D, Biswal DK. *Complete plastid ndh gene retention in* Nepenthes khasiana*: a comparative plastome analysis across carnivorous plant lineages.* Bioinformatics Advances (in revision, BIOADV-2026-545).

The study compares 49 carnivorous plant plastomes from nine families and five orders, with *Amborella trichopoda* as the only non-carnivorous taxon (outgroup). The focal plastomes are *Nepenthes khasiana* (MH923233 / NC_051455) and *N. gracilis* (PP035775).

## Contents

| Path | What it is |
| --- | --- |
| `data/plastome_accessions.csv` | The 50 plastomes analysed: species, family, role, GenBank accession, length (Supplementary Table S1) |
| `data/ndh_status_table1.csv` | ndhA–ndhK status (F functional, P pseudogene, – absent) for the Table 1 lineages |
| `scripts/fetch_plastomes.py` | Downloads the 50 GenBank records from NCBI |
| `scripts/ndh_status.py` | Scores ndhA–ndhK in every plastome from its annotation and assigns ndh-complete / ndh-partial / ndh-null classes |
| `scripts/check_ndh_matrix.py` | Checks that every row of a curated ndh matrix accounts for all 11 loci |
| `scripts/rscu_pca.py` | RSCU from concatenated CDS, then PCA on the centred and scaled matrix (cross-check of Figure 5) |
| `scripts/kmer_distance.py` | 15-mer Mash-style distances and a UPGMA tree (cross-check of Figure 3) |
| `scripts/genbank_utils.py` | Minimal GenBank reader used by the scripts (no Biopython needed) |
| `scripts/original/` | The authors' original MEGA, DnaSP, MAUVE, PGR-TK and PGGB files and commands |

## Requirements

Python ≥ 3.9 with `numpy`, `scipy` and `pandas`.

```bash
pip install numpy scipy pandas
```

The phylogeny (MEGA X), RSCU (MEGA 12.1), dN/dS (DnaSP v6), synteny (progressiveMauve) and pangenome graphs (PGR-TK, PGGB) were run with those programs; their settings and outputs are in `scripts/original/`.

## Reproducing the comparative analyses

```bash
# 1. download the plastomes (about 8 MB)
python scripts/fetch_plastomes.py --email you@example.org

# 2. ndh status of all 50 plastomes, from the GenBank annotations
python scripts/ndh_status.py data/genbank/*.gb --out results/ndh_status.csv

# 3. check a curated matrix (every row must sum to 11 loci)
python scripts/check_ndh_matrix.py data/ndh_status_table1.csv

# 4. codon usage and PCA on the 49 carnivorous plastomes
python scripts/rscu_pca.py data/genbank/*.gb --outdir results/codon            # 59 codons
python scripts/rscu_pca.py data/genbank/*.gb --outdir results/codon --keep-trp # 60 codons

# 5. alignment-free 15-mer clustering
python scripts/kmer_distance.py data/genbank/*.gb --k 15 --outdir results/kmer
```

`ndh_status.py` scores loci automatically from the annotation (open reading frame present, no internal stop, at least 80% of the reference length, not flagged `/pseudo`). The published matrix (Figure 3, Table 1) was curated by hand, so treat any differences as cases to inspect, not as corrections.

## ndh retention classes

| Class | Definition |
| --- | --- |
| ndh-complete | All 11 plastid ndh loci functional |
| ndh-partial | 1–10 loci functional; the rest pseudogenized or absent |
| ndh-null | No functional locus |

These classes describe plastome gene content only. They are not measures of photosynthetic or nutritional physiology.

## Licence

MIT (code). The GenBank records belong to their submitters.
