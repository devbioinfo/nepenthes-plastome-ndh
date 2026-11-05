# nepenthes-plastome-ndh

Data, code and results for:

> Konhar R, Dash D, Biswal DK. *A comparative pan-plastome analysis of carnivorous plants reveals structural conservatism in* Nepenthes khasiana *and lineage-specific loss and pseudogenization of* ndh *genes.* Bioinformatics Advances (in revision, BIOADV-2026-545).

The dataset is 49 carnivorous plant plastomes from nine families and five orders, plus *Amborella trichopoda* as the only non-carnivorous taxon (outgroup). The focal plastomes are *Nepenthes khasiana* (MH923233 / NC_051455) and *N. gracilis* (PP035775).

## Reproduce everything

```bash
bash run_all.sh you@example.org     # e-mail is passed to NCBI E-utilities, as NCBI requests
```

Requirements: Python ≥ 3.9 with numpy, scipy, pandas, matplotlib, biopython and openpyxl; MAFFT v7.505; IQ-TREE v2.0.7; BLAST+ v2.12.0. On Debian/Ubuntu, `apt-get install mafft iqtree ncbi-blast+` followed by `pip install numpy scipy pandas matplotlib biopython openpyxl` is enough. The run takes about 20 minutes on two cores, most of it IQ-TREE.

## What each step does

| Step | Script | Output | Manuscript |
| --- | --- | --- | --- |
| Download the 50 GenBank records | `scripts/fetch_plastomes.py` | `data/genbank/` | Table S1 |
| Score *ndhA–ndhK* from annotations (intact ORF ≥ 80% of reference; ACG/GTG start extension) | `scripts/ndh_status.py` | `results/ndh_status_auto*.csv` | Methods 2.5 |
| Search unannotated loci (tblastn, ≥ 50% identity) and test ORFs | `scripts/ndh_tblastn_check.py`, `scripts/ndh_resolve_unannotated.py` | `results/ndh_tblastn_absent_loci.csv`, `results/ndh_unannotated_resolved.csv`, `results/ndh_pseudogene_recheck.csv` | Methods 2.5 |
| Curated ndh matrix and classes | `scripts/ndh_curate.py`, `scripts/check_ndh_matrix.py` | `results/ndh_status_curated.csv` | Figure 3, Table 1, Table S2 |
| Codon alignment of 68 shared protein-coding genes (MAFFT L-INS-i) | `scripts/build_cds_alignment.py` | `results/phylogeny/concatenated_cds.fasta` (56,523 sites) | Methods 2.3, Table S3 |
| ML tree (ModelFinder → GTR+F+I+G4 by AICc; 1,000 UFBoot + SH-aLRT) | IQ-TREE call in `run_all.sh` | `results/phylogeny/cds_ml.*` | Figure 2A |
| 15-mer Mash distances and neighbour-joining tree | `scripts/kmer_distance.py` | `results/kmer/` | Figure 3 |
| Whole-plastome ANIb (1,020-nt fragments, dc-megablast) and UPGMA | `scripts/plastome_ani.py` | `results/ani/` | Figure 6 |
| RSCU (59 codons) and PCA | `scripts/rscu_pca.py` | `results/codon/` | Figure 5 |
| ENC and GC3s vs. Wright (1990) expectation | `scripts/enc_gc3.py` | `results/codon/enc_gc3.csv` | Figure 5C |
| dN, dS, dN/dS vs *N. khasiana* (Nei–Gojobori, Biopython) | `scripts/dnds_ng86.py` | `results/dnds/` | Figure S1 |
| Figures | `scripts/figures.py`, `scripts/compose_fig2.py` | `figures/` | Figures 2, 3, 5, 6, S1 |

Figure 2 panels B and C are taken unchanged from the original figure; `compose_fig2.py` adds them to the new panel A. Figures 1 (OGDRAW maps) and 4 (MAUVE/pangenome graph) were not regenerated; their settings are given in Methods 2.2 and 2.4.

## ndh retention classes

| Class | Definition | Plastomes |
| --- | --- | --- |
| ndh-complete | all 11 loci functional | 21 carnivorous + *Amborella* |
| ndh-partial | 1–10 loci functional | 9 |
| ndh-null | no functional locus | 19 |

The classes describe plastome gene content only, not photosynthetic or nutritional physiology. `results/ndh_curation_changes.csv` lists every locus where sequence evidence changed the annotation-based call. Examples are *Triphyophyllum ndhB* and *Utricularia macrorhiza ndhA/F/G/K*, which are intact but not annotated in GenBank.

## Licence

MIT (code). GenBank records belong to their submitters.
