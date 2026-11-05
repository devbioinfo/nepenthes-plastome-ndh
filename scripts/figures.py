#!/usr/bin/env python3
"""Draw Figures 2A, 3, 5, 6 and Supplementary Figure S1 from the results/ tables.

Usage: python scripts/figures.py --outdir figures
"""
import argparse, os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import squareform
from Bio import Phylo
from io import StringIO

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7, 'axes.linewidth': 0.6,
                     'axes.edgecolor': '#52514e', 'xtick.color': '#52514e', 'ytick.color': '#52514e',
                     'axes.labelcolor': '#0b0b0b', 'text.color': '#0b0b0b', 'savefig.dpi': 300})
INK, INK2, GRID = '#0b0b0b', '#52514e', '#d9d8d3'
ORDER_COL = {'Caryophyllales': '#2a78d6', 'Lamiales': '#eb6834', 'Ericales': '#1baf7a',
             'Oxalidales': '#eda100', 'Poales': '#e87ba4', 'Amborellales': '#8a8984'}
FAM_ORDER = {'Droseraceae': 'Caryophyllales', 'Drosophyllaceae': 'Caryophyllales', 'Dioncophyllaceae': 'Caryophyllales',
             'Nepenthaceae': 'Caryophyllales', 'Lentibulariaceae': 'Lamiales', 'Sarraceniaceae': 'Ericales',
             'Roridulaceae': 'Ericales', 'Cephalotaceae': 'Oxalidales', 'Bromeliaceae': 'Poales', 'Amborellaceae': 'Amborellales'}
FAM_MARK = {'Droseraceae': 'o', 'Drosophyllaceae': 's', 'Dioncophyllaceae': '^', 'Nepenthaceae': 'D',
            'Lentibulariaceae': 'o', 'Sarraceniaceae': 'o', 'Roridulaceae': 's', 'Cephalotaceae': 'o',
            'Bromeliaceae': 'o', 'Amborellaceae': 'o'}
TRAP = {'Aldrovanda': 'snap', 'Dionaea': 'snap', 'Drosera': 'adhesive', 'Drosophyllum': 'adhesive',
        'Triphyophyllum': 'adhesive', 'Nepenthes': 'pitfall', 'Sarracenia': 'pitfall', 'Heliamphora': 'pitfall',
        'Darlingtonia': 'pitfall', 'Cephalotus': 'pitfall', 'Brocchinia': 'pitfall', 'Roridula': 'adhesive',
        'Pinguicula': 'adhesive', 'Genlisea': 'eel', 'Utricularia': 'suction', 'Amborella': 'outgroup'}
NDH = ['ndhA', 'ndhB', 'ndhC', 'ndhD', 'ndhE', 'ndhF', 'ndhG', 'ndhH', 'ndhI', 'ndhJ', 'ndhK']


def meta():
    acc = pd.read_csv('data/plastome_accessions.csv')
    acc['taxon'] = acc.species.str.replace(' x ', '_x_').str.replace(' ', '_')
    acc['order'] = acc.family.map(FAM_ORDER)
    acc['trap'] = acc.species.str.split().str[0].map(TRAP)
    return acc.set_index('taxon')


def pretty(t):
    if t.startswith('Nepenthes_ventricosa'):
        return 'Nepenthes ventricosa × N. alata'
    return t.replace('_x_', ' × ').replace('_', ' ')


def tree_layout(newick, outgroup='Amborella_trichopoda'):
    tr = Phylo.read(StringIO(open(newick).read()), 'newick')
    tr.root_with_outgroup(outgroup)
    tr.ladderize(reverse=True)
    terms = tr.get_terminals()
    y = {c: i for i, c in enumerate(terms)}
    depths = tr.depths() if any(c.branch_length for c in tr.find_clades()) else tr.depths(unit_branch_lengths=True)

    def sety(c):
        if c in y:
            return y[c]
        ys = [sety(ch) for ch in c.clades]
        y[c] = (min(ys) + max(ys)) / 2
        return y[c]
    sety(tr.root)
    return tr, terms, y, depths


def draw_tree(ax, tr, y, depths, support=False, lw=0.7, scale_label=None):
    for c in tr.find_clades():
        x1 = depths[c]
        for ch in c.clades:
            ax.plot([x1, depths[ch]], [y[ch], y[ch]], color=INK2, lw=lw, solid_capstyle='butt')
        if c.clades:
            ys = [y[ch] for ch in c.clades]
            ax.plot([x1, x1], [min(ys), max(ys)], color=INK2, lw=lw)
            if support and c.confidence is not None and c is not tr.root:
                pass
    ax.invert_yaxis()
    ax.axis('off')


# ---------------- Figure 3 ----------------
def fig3(out):
    m = meta()
    tr, terms, y, depths = tree_layout('results/kmer/kmer_nj.nwk')
    ndh = pd.read_csv('results/ndh_status_curated.csv')
    ndh['taxon'] = ndh.species.str.replace(' x ', '_x_').str.replace(' ', '_')
    ndh = ndh.set_index('taxon')
    n = len(terms)
    fig = plt.figure(figsize=(7.2, 7.6))
    axT = fig.add_axes([0.02, 0.06, 0.30, 0.86])
    axL = fig.add_axes([0.32, 0.06, 0.25, 0.86], sharey=axT)
    axM = fig.add_axes([0.57, 0.06, 0.27, 0.86], sharey=axT)
    axS = fig.add_axes([0.85, 0.06, 0.12, 0.86], sharey=axT)
    draw_tree(axT, tr, y, depths)
    xmax = max(depths[t] for t in terms)
    for t in terms:
        axT.plot([depths[t], xmax], [y[t], y[t]], color=GRID, lw=0.4, ls=(0, (1, 2)))
    axT.set_xlim(-0.01, xmax * 1.01)
    axT.plot([0, 0.05], [n + 0.2, n + 0.2], color=INK, lw=0.8, clip_on=False)
    axT.text(0.025, n + 1.4, '0.05', ha='center', va='center', fontsize=6, color=INK2, clip_on=False)
    for t in terms:
        o = m.loc[t.name, 'order']
        axL.add_patch(Rectangle((0.0, y[t] - 0.35), 0.035, 0.7, color=ORDER_COL[o], lw=0))
        st = 'italic'
        axL.text(0.06, y[t], pretty(t.name), va='center', fontsize=6.2, style=st, color=INK)
    axL.axis('off'); axL.set_xlim(0, 1)
    sym = {'F': dict(fc=INK, ec=INK), 'P': dict(fc='white', ec=INK), '-': dict(fc='#f0efec', ec='#f0efec')}
    for t in terms:
        r = ndh.loc[t.name]
        for j, g in enumerate(NDH):
            axM.add_patch(Rectangle((j + 0.12, y[t] - 0.38), 0.76, 0.76, lw=0.6, **sym[r[g]]))
    for j, g in enumerate(NDH):
        axM.text(j + 0.5, -1.0, g, rotation=90, ha='center', va='bottom', fontsize=6, style='italic', color=INK)
    axM.set_xlim(0, 11); axM.axis('off')
    for t in terms:
        L = m.loc[t.name, 'length_bp'] / 1000
        axS.barh(y[t], L - 100, left=100, height=0.6, color='#86b6ef', lw=0)
    axS.set_xlim(100, 200); axS.set_xticks([100, 150, 200])
    axS.tick_params(axis='y', left=False, labelleft=False); axS.tick_params(axis='x', labelsize=6, length=2)
    for s in ('top', 'right', 'left'):
        axS.spines[s].set_visible(False)
    axS.set_xlabel('Plastome size (kb)', fontsize=6.5)
    axS.text(100, -1.0, ' ', fontsize=6)
    handles = [Patch(fc=INK, ec=INK, label='functional'), Patch(fc='white', ec=INK, label='pseudogene'),
               Patch(fc='#f0efec', ec='#f0efec', label='absent')]
    handles += [Patch(fc=c, label=o) for o, c in ORDER_COL.items()]
    fig.legend(handles=handles, loc='upper center', ncol=9, frameon=False, fontsize=6, bbox_to_anchor=(0.5, -0.005),
               handlelength=1.0, columnspacing=1.0)
    fig.text(0.02, 0.965, '15-mer neighbour-joining tree', fontsize=7, color=INK2)
    fig.text(0.57, 0.965, 'Plastid ndh genes', fontsize=7, color=INK2)
    fig.savefig(out, bbox_inches='tight'); plt.close(fig)
    return [t.name for t in terms]


# ---------------- Figure 5 ----------------
def fig5(out):
    m = meta()
    r = pd.read_csv('results/codon/rscu.csv', index_col=0)
    r.index = r.index.str.replace(' x ', '_x_').str.replace(' ', '_')
    sc = pd.read_csv('results/codon/pca_scores.csv', index_col=0)
    sc.index = sc.index.str.replace(' x ', '_x_').str.replace(' ', '_')
    var = pd.read_csv('results/codon/pca_variance.csv', index_col=0).iloc[:, 0]
    famorder = ['Droseraceae', 'Drosophyllaceae', 'Dioncophyllaceae', 'Nepenthaceae', 'Lentibulariaceae',
                'Sarraceniaceae', 'Roridulaceae', 'Cephalotaceae', 'Bromeliaceae']
    rows = sorted(r.index, key=lambda t: (famorder.index(m.loc[t, 'family']), t))
    B = 'TCAG'; AA = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    code = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate(B) for j, b in enumerate(B) for k, c in enumerate(B)}
    cols = sorted(r.columns, key=lambda c: (code[c], c))
    M = r.loc[rows, cols].values
    fig = plt.figure(figsize=(7.2, 4.6))
    ax = fig.add_axes([0.13, 0.12, 0.47, 0.82])
    cmap = LinearSegmentedColormap.from_list('div', ['#2a78d6', '#f0efec', '#eb6834'])
    im = ax.imshow(M, aspect='auto', cmap=cmap, norm=TwoSlopeNorm(vmin=0, vcenter=1, vmax=max(2.5, np.nanmax(M))),
                   interpolation='nearest')
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([pretty(t) for t in rows], fontsize=4.6, style='italic')
    ax.set_xticks(range(len(cols))); ax.set_xticklabels([f'{c}' for c in cols], rotation=90, fontsize=4.2, family='DejaVu Sans Mono')
    for k, c in enumerate(cols):
        if k == 0 or code[cols[k - 1]] != code[c]:
            ax.text(k, -1.2, code[c], ha='left', va='bottom', fontsize=5, color=INK2)
    ax.tick_params(length=0)
    cax = fig.add_axes([0.13, 0.035, 0.2, 0.015])
    cb = fig.colorbar(im, cax=cax, orientation='horizontal'); cb.ax.tick_params(labelsize=5.5, length=2)
    cb.set_label('RSCU (1 = no bias)', fontsize=6)
    fig.text(0.13, 0.975, 'A', fontsize=10, weight='bold')
    ax2 = fig.add_axes([0.68, 0.32, 0.30, 0.58])
    for t in sc.index:
        f = m.loc[t, 'family']
        ax2.scatter(sc.loc[t, 'PC1'], sc.loc[t, 'PC2'], s=22, marker=FAM_MARK[f], color=ORDER_COL[FAM_ORDER[f]],
                    edgecolor='white', linewidth=0.6, zorder=3)
    for t in ['Nepenthes_khasiana']:
        ax2.annotate('N. khasiana', (sc.loc[t, 'PC1'], sc.loc[t, 'PC2']), xytext=(8, -12), textcoords='offset points',
                     fontsize=6, style='italic', arrowprops=dict(arrowstyle='-', lw=0.5, color=INK2))
    ax2.axhline(0, color=GRID, lw=0.5, zorder=0); ax2.axvline(0, color=GRID, lw=0.5, zorder=0)
    ax2.set_xlabel(f'PC1 ({var.iloc[0]:.1%})'); ax2.set_ylabel(f'PC2 ({var.iloc[1]:.1%})')
    for s in ('top', 'right'):
        ax2.spines[s].set_visible(False)
    ax2.tick_params(labelsize=6, length=2)
    hs = [Line2D([], [], ls='', marker=FAM_MARK[f], color=ORDER_COL[FAM_ORDER[f]], markeredgecolor='white',
                 markersize=5, label=f) for f in famorder]
    fig.legend(handles=hs, loc='lower left', bbox_to_anchor=(0.66, 0.02), ncol=2, frameon=False, fontsize=6,
               handletextpad=0.3, columnspacing=0.8)
    fig.text(0.63, 0.975, 'B', fontsize=10, weight='bold')
    fig.savefig(out, bbox_inches='tight'); plt.close(fig)


# ---------------- Figure 6 ----------------
def fig6(out):
    A = pd.read_csv('results/ani/ani_matrix.tsv', sep='\t', index_col=0)
    D = 100 - A.values
    np.fill_diagonal(D, 0); D = (D + D.T) / 2
    Z = linkage(squareform(D, checks=False), method='average')
    order = leaves_list(Z)
    names = A.index[order]
    M = A.values[np.ix_(order, order)]
    fig = plt.figure(figsize=(7.2, 7.0))
    axd = fig.add_axes([0.20, 0.80, 0.62, 0.12])
    dendrogram(Z, ax=axd, color_threshold=0, above_threshold_color=INK2, no_labels=True)
    axd.axis('off')
    ax = fig.add_axes([0.20, 0.06, 0.62, 0.73])
    cmap = LinearSegmentedColormap.from_list('seq', ['#cde2fb', '#86b6ef', '#3987e5', '#1c5cab', '#0d366b'])
    im = ax.imshow(M, cmap=cmap, vmin=80, vmax=100, interpolation='nearest', aspect='auto')
    ax.set_xticks(range(len(names))); ax.set_xticklabels([pretty(t) for t in names], rotation=90, fontsize=4.6, style='italic')
    ax.set_yticks(range(len(names))); ax.set_yticklabels([pretty(t) for t in names], fontsize=4.6, style='italic')
    ax.tick_params(length=0)
    cax = fig.add_axes([0.85, 0.30, 0.015, 0.3])
    cb = fig.colorbar(im, cax=cax); cb.ax.tick_params(labelsize=6, length=2)
    cb.set_label('Average nucleotide identity (%)', fontsize=6.5)
    fig.savefig(out, bbox_inches='tight'); plt.close(fig)


# ---------------- Figure S1 ----------------
def figS1(out):
    m = meta()
    d = pd.read_csv('results/dnds/dnds_vs_khasiana.csv')
    d['family'] = d.species.map(m.family)
    d = d[d.flag.isna()]
    d['ratio'] = pd.to_numeric(d.dN_dS, errors='coerce')
    fams = ['Nepenthaceae', 'Dioncophyllaceae', 'Drosophyllaceae', 'Droseraceae', 'Sarraceniaceae', 'Roridulaceae',
            'Cephalotaceae', 'Lentibulariaceae', 'Bromeliaceae', 'Amborellaceae']
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.0), sharey=True)
    colc = {'core': '#2a78d6', 'ndh': '#eb6834'}
    rng = np.random.default_rng(1)
    for ax, (col, lab) in zip(axes, [('dN', 'dN'), ('dS', 'dS'), ('ratio', 'dN/dS')]):
        for i, f in enumerate(fams):
            for k, cl in enumerate(['core', 'ndh']):
                v = d[(d.family == f) & (d['class'] == cl)][col].dropna().values
                if not len(v):
                    continue
                yy = i + (-0.17 if cl == 'core' else 0.17)
                ax.scatter(v, yy + rng.uniform(-0.09, 0.09, len(v)), s=5, color=colc[cl], alpha=0.55, lw=0)
                ax.plot([np.median(v)] * 2, [yy - 0.16, yy + 0.16], color=INK, lw=1.2)
        ax.set_xlabel(lab, style='italic' if lab != 'dN/dS' else 'normal')
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
        ax.tick_params(labelsize=6, length=2)
        ax.grid(axis='x', color=GRID, lw=0.4); ax.set_axisbelow(True)
    axes[0].set_yticks(range(len(fams))); axes[0].set_yticklabels(fams, fontsize=6.5)
    axes[0].invert_yaxis()
    axes[2].set_xlim(0, 1.0)
    hs = [Line2D([], [], ls='', marker='o', color=colc['core'], label='core photosynthetic genes (11)'),
          Line2D([], [], ls='', marker='o', color=colc['ndh'], label='plastid ndh genes (11)'),
          Line2D([], [], color=INK, lw=1.2, label='median')]
    fig.legend(handles=hs, loc='lower center', ncol=3, frameon=False, fontsize=6.5, bbox_to_anchor=(0.55, -0.04))
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(out, bbox_inches='tight'); plt.close(fig)


def fig2a(out):
    m = meta()
    tr, terms, y, depths = tree_layout('results/phylogeny/cds_ml.treefile')
    fig = plt.figure(figsize=(4.6, 7.2))
    ax = fig.add_axes([0.02, 0.04, 0.52, 0.92])
    draw_tree(ax, tr, y, depths, lw=0.8)
    xmax = max(depths[t] for t in terms)
    for c in tr.get_nonterminals():
        lab = c.name or ''
        if c is tr.root or '/' not in lab:
            continue
        ufb = float(lab.split('/')[1])
        if ufb < 100:
            ax.text(depths[c] - xmax * 0.005, y[c] - 0.35, f'{int(round(ufb))}', ha='right', va='bottom',
                    fontsize=5, color=INK2)
    for t in terms:
        o = m.loc[t.name, 'order']
        ax.text(depths[t] + xmax * 0.012, y[t], pretty(t.name), va='center', fontsize=5.6, style='italic')
        ax.add_patch(Rectangle((xmax * 1.50, y[t] - 0.4), xmax * 0.03, 0.8, color=ORDER_COL[o], lw=0, clip_on=False))
        ax.text(xmax * 1.55, y[t], m.loc[t.name, 'trap'], va='center', fontsize=5.4, color=INK2, clip_on=False)
    ax.set_xlim(0, xmax * 1.49)
    ax.plot([0, 0.05], [len(terms) + 0.5] * 2, color=INK, lw=0.8, clip_on=False)
    ax.text(0.025, len(terms) + 1.6, '0.05 substitutions/site', ha='left', fontsize=5.5, color=INK2, clip_on=False)
    # family brackets
    fams = {}
    for t in terms:
        fams.setdefault(m.loc[t.name, 'family'], []).append(y[t])
    for f, ys in fams.items():
        ax.plot([xmax * 1.78] * 2, [min(ys) - 0.3, max(ys) + 0.3], color=INK2, lw=0.6, clip_on=False)
        ax.text(xmax * 1.81, np.mean(ys), f, va='center', fontsize=5.8, clip_on=False)
    hs = [Patch(fc=c, label=o) for o, c in ORDER_COL.items()]
    fig.legend(handles=hs, loc='lower center', ncol=3, frameon=False, fontsize=5.8, bbox_to_anchor=(0.5, -0.02))
    fig.savefig(out, bbox_inches='tight'); plt.close(fig)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--outdir', default='figures'); ap.add_argument('--only', default='')
    a = ap.parse_args(); os.makedirs(a.outdir, exist_ok=True)
    todo = a.only.split(',') if a.only else ['3', '5', '6', 'S1', '2A']
    if '3' in todo: fig3(os.path.join(a.outdir, 'Figure3.png'))
    if '5' in todo: fig5(os.path.join(a.outdir, 'Figure5.png'))
    if '6' in todo: fig6(os.path.join(a.outdir, 'Figure6.png'))
    if 'S1' in todo: figS1(os.path.join(a.outdir, 'FigureS1.png'))
    if '2A' in todo and os.path.exists('results/phylogeny/cds_ml.contree'): fig2a(os.path.join(a.outdir, 'Figure2A.png'))
