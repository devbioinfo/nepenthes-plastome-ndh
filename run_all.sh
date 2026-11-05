#!/usr/bin/env bash
# Reproduce every comparative analysis and figure from the GenBank accessions.
# Requirements: python>=3.9 (numpy, scipy, pandas, matplotlib, biopython, openpyxl),
#               MAFFT v7.505, IQ-TREE v2.0.7, BLAST+ v2.12.0.
set -euo pipefail
EMAIL=${1:?usage: bash run_all.sh your@email.org}
python3 scripts/fetch_plastomes.py --email "$EMAIL"
GB=(data/genbank/*.gb)

# 1. ndh status (Methods 2.5)
python3 scripts/ndh_status.py "${GB[@]}" --out results/ndh_status_auto.csv
python3 scripts/ndh_tblastn_check.py results/ndh_status_auto.csv data/genbank > results/ndh_tblastn_absent_loci.csv
python3 scripts/ndh_resolve_unannotated.py results/ndh_tblastn_absent_loci.csv results/ndh_status_auto.csv data/genbank > results/ndh_unannotated_resolved.csv
python3 - <<'PY'
import csv
NDH=['ndhA','ndhB','ndhC','ndhD','ndhE','ndhF','ndhG','ndhH','ndhI','ndhJ','ndhK']
with open('results/_pseudo_loci.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['species','gene','query_coverage','call'])
    for r in csv.DictReader(open('results/ndh_status_auto.csv')):
        for g in NDH:
            if r[g]=='P': w.writerow([r['species'],g,'','CHECK'])
PY
python3 scripts/ndh_resolve_unannotated.py results/_pseudo_loci.csv results/ndh_status_auto.csv data/genbank > results/ndh_pseudogene_recheck.csv
python3 scripts/ndh_curate.py
python3 scripts/check_ndh_matrix.py results/ndh_status_curated.csv > results/ndh_check.txt

# 2. ML phylogeny (Methods 2.3)
python3 scripts/build_cds_alignment.py "${GB[@]}" --outdir results/phylogeny
(cd results/phylogeny && iqtree2 -s concatenated_cds.fasta -m MFP -mset JC,F81,K2P,HKY,SYM,GTR \
   -mrate E,I,G,I+G -merit AICc -B 1000 -alrt 1000 -o Amborella_trichopoda -T AUTO -seed 12345 --prefix cds_ml -redo)

# 3. k-mer tree, ANI, codon usage, dN/dS (Methods 2.3, 2.6, 2.7)
python3 scripts/kmer_distance.py "${GB[@]}" --k 15 --outdir results/kmer
python3 scripts/plastome_ani.py "${GB[@]}" --outdir results/ani
python3 scripts/rscu_pca.py "${GB[@]}" --outdir results/codon
python3 scripts/enc_gc3.py "${GB[@]}" --out results/codon/enc_gc3.csv
python3 scripts/dnds_ng86.py "${GB[@]}" --curated results/ndh_status_curated.csv --outdir results/dnds

# 4. figures
python3 scripts/figures.py --outdir figures
echo "Figure 2 panels B and C are taken from the original figure: python3 scripts/compose_fig2.py <original_figure2.png>"
