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


# ---------------- Figure 2A, classic style (collapsed clades, original colour coding and trap icons) ----------------
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
ICON = {k: plt.imread(f'figures/icons/{k}.png') for k in ('eel', 'suction', 'flypaper', 'pitfall', 'snap')}
TRAPNAME = {'eel': 'Eel (lobster-pot) trap', 'suction': 'Suction (bladder) trap', 'flypaper': 'Adhesive (flypaper) trap',
            'pitfall': 'Pitfall trap', 'tank': 'Pitfall (tank) trap', 'snap': 'Snap trap'}
ICON['tank'] = ICON['pitfall']; ICON['snap_d'] = ICON['snap']; ICON['fly_d'] = ICON['flypaper']
TRAPNAME['snap_d'] = r'Snap trap ($\it{Dionaea}$, $\it{Aldrovanda}$)'
TRAPNAME['fly_d'] = r'Adhesive (flypaper) trap ($\it{Drosera}$)'
# (family label, genus prefixes, triangle fill, text colour, traps) - colours follow the original figure
CLADES = [
    ('Lentibulariaceae', ('Utricularia',), '#ddd0d0', '#1a1a1a', ['suction']),
    ('Lentibulariaceae', ('Genlisea',), '#d6f3fb', '#1f9fd0', ['eel']),
    ('Lentibulariaceae', ('Pinguicula',), '#d7f0c6', '#5d9e3a', ['flypaper']),
    ('Sarraceniaceae', ('Sarracenia', 'Heliamphora', 'Darlingtonia'), '#d9d1f4', '#5b4bc4', ['pitfall']),
    ('Droseraceae', ('Drosera', 'Dionaea', 'Aldrovanda'), '#f9d3d9', '#b03a4e', ['snap_d', 'fly_d']),
    ('Nepenthaceae', ('Nepenthes',), '#f7aa2c', '#e07f00', ['pitfall']),
    ('Bromeliaceae', ('Brocchinia',), '#ea4a1e', '#c8300c', ['tank']),
]
SINGLE = {'Roridula_gorgonias': ('Roridulaceae', '#a926b0', ['flypaper']),
          'Drosophyllum_lusitanicum': ('Drosophyllaceae', '#1f4e9e', ['flypaper']),
          'Triphyophyllum_peltatum': ('Dioncophyllaceae', '#1f4e9e', ['flypaper']),
          'Cephalotus_follicularis': ('Cephalotaceae', '#9a4f1c', ['pitfall']),
          'Amborella_trichopoda': ('Amborellaceae (outgroup)', '#2e7d32', [])}


def fig2a_classic(out):
    tr = Phylo.read('results/phylogeny/cds_ml.treefile', 'newick')
    tr.root_with_outgroup('Amborella_trichopoda'); tr.ladderize(reverse=True)
    depths = tr.depths()
    parent = {id(ch): c for c in tr.find_clades() for ch in c.clades}
    mrca = {}
    for lab, pre, fill, tcol, traps in CLADES:
        tips = [x for x in tr.get_terminals() if x.name.startswith(pre)]
        mrca[id(tr.common_ancestor(tips))] = (lab, fill, tcol, traps)
    y = {}; rows = []; cur = [0.0]

    def place(c):
        if id(c) in mrca:
            order = c.get_terminals()
            h = max(len(order), 2) * 0.62 + 0.6
            y[c] = cur[0] + h / 2; rows.append(('clade', c, cur[0], cur[0] + h, order)); cur[0] += h + 0.6
            return y[c]
        if c.is_terminal():
            y[c] = cur[0] + 0.6; rows.append(('tip', c, cur[0], cur[0] + 1.2, [c])); cur[0] += 1.9
            return y[c]
        ys = [place(ch) for ch in c.clades]
        y[c] = (min(ys) + max(ys)) / 2
        return y[c]
    place(tr.root)
    H = cur[0]; xmax = max(depths.values())
    fig = plt.figure(figsize=(7.4, 7.8))
    ax = fig.add_axes([0.02, 0.05, 0.96, 0.92])
    stemcol = {}
    for kind, c, *_ in rows:
        stemcol[id(c)] = (mrca[id(c)][2] if kind == 'clade' else SINGLE[c.name][1])

    def sup(c):
        lab = c.name or ''
        if '/' not in lab or c is tr.root:
            return
        p = parent.get(id(c))
        x = (depths[p] + depths[c]) / 2 if p is not None else depths[c]
        ax.text(x, y[c] - 0.22, f"{float(lab.split('/')[1]):.0f}%", ha='center', va='bottom', fontsize=4.8, color=INK2,
                bbox=dict(boxstyle='square,pad=0.05', fc='white', ec='none'), zorder=4)

    def draw(c):
        if id(c) in mrca or c.is_terminal():
            return
        x1 = depths[c]; ys = [y[ch] for ch in c.clades]
        ax.plot([x1, x1], [min(ys), max(ys)], color=INK, lw=0.9)
        for ch in c.clades:
            col = stemcol.get(id(ch), INK)
            ax.plot([x1, depths[ch]], [y[ch], y[ch]], color=col if col != '#1a1a1a' else INK, lw=1.1 if id(ch) in stemcol else 0.9)
            draw(ch)
        sup(c); ax.plot(x1, y[c], 'o', ms=1.8, color=INK)
    draw(tr.root)

    def icon(name, xy, zoom=0.2):
        ax.add_artist(AnnotationBbox(OffsetImage(ICON[name], zoom=zoom), xy, frameon=False, box_alignment=(0, 0.5)))

    lx = xmax * 1.04; ix = lx + xmax * 0.64; fx = ix + xmax * 0.74
    def traps_at(traps, yc):
        n = len(traps)
        for k, tp in enumerate(traps):
            yy = yc + (k - (n - 1) / 2) * 1.35
            icon(tp, (ix, yy)); ax.text(ix + xmax * 0.10, yy, TRAPNAME[tp], va='center', fontsize=5.8, color='#1c3f94')
    for kind, c, y0, y1, tips in rows:
        if kind == 'clade':
            lab, fill, tcol, traps = mrca[id(c)]
            x0 = depths[c]; xe = max(depths[t] for t in tips)
            ax.add_patch(plt.Polygon([[x0, y[c]], [xe, y0 + 0.2], [xe, y1 - 0.2]], closed=True, fc=fill, ec=fill, lw=0.6))
            sup(c); ax.plot(x0, y[c], 'o', ms=1.8, color=INK)
            step = (y1 - y0 - 0.6) / max(len(tips) - 1, 1)
            for k, t in enumerate(tips):
                ax.text(lx, y0 + 0.3 + k * step, pretty(t.name), va='center', fontsize=5.6, style='italic', color=tcol)
            ax.plot([fx - xmax * 0.03] * 2, [y0 + 0.2, y1 - 0.2], color=tcol, lw=0.8)
            ax.text(fx, y[c], lab, va='center', fontsize=7, weight='bold', color=tcol)
            traps_at(traps, y[c])
        else:
            fam, col, traps = SINGLE[c.name]
            ax.plot(depths[c], y[c], 's', ms=3.2, color=col)
            ax.text(lx, y[c], pretty(c.name), va='center', fontsize=5.6, style='italic', color=col)
            ax.text(fx, y[c], fam, va='center', fontsize=7, weight='bold', color=col)
            traps_at(traps, y[c])
    for kind, c, y0, y1, tips in rows:
        xe = max(depths[t] for t in tips)
        ax.plot([xe + xmax * 0.01, lx - xmax * 0.01], [y[c], y[c]], color=GRID, lw=0.5, ls=(0, (1, 2)))
    ax.set_xlim(-xmax * 0.02, xmax * 2.95); ax.set_ylim(H + 1.5, -0.8)
    ax.plot([0, 0.05], [H + 0.6] * 2, color=INK, lw=0.8)
    ax.text(0.0, H + 1.4, '0.05 substitutions/site', fontsize=5.5, color=INK2)
    ax.axis('off')
    fig.savefig(out, bbox_inches='tight'); plt.close(fig)


# ---------------- Figure 3, illustrated style (cladogram, clade bands, plant drawings) ----------------
from matplotlib.colors import to_rgb

FAM3 = [  # (key, genus prefixes, branch/label colour, band colour or None, family label, label colour)
    ('Genlisea', ('Genlisea',), '#1f6fe0', None, None, None),
    ('Pinguicula', ('Pinguicula',), '#d60000', None, None, None),
    ('Utricularia', ('Utricularia',), '#157a15', None, None, None),
    ('Lentibulariaceae', ('Genlisea', 'Pinguicula', 'Utricularia'), None, '#7fd46a', 'Lentibulariaceae', '#2e9e45'),
    ('Sarraceniaceae', ('Sarracenia', 'Heliamphora', 'Darlingtonia'), '#4646d8', None, 'Sarraceniaceae', '#4a6fb5'),
    ('Roridula', ('Roridula',), '#e0128c', None, 'Roridulaceae', '#e0128c'),
    ('Ericales', ('Sarracenia', 'Heliamphora', 'Darlingtonia', 'Roridula'), None, '#e8418f', None, None),
    ('Droseraceae', ('Drosera', 'Dionaea', 'Aldrovanda'), '#1aa6bd', '#5fe0ee', 'Droseraceae', '#2bbfd6'),
    ('Nepenthaceae', ('Nepenthes',), '#7a0b0b', '#cf6a32', 'Nepenthaceae', '#7a0b0b'),
    ('Drosophyllaceae', ('Drosophyllum',), '#1030e8', None, 'Drosophyllaceae', '#1030e8'),
    ('Dioncophyllaceae', ('Triphyophyllum',), '#1030e8', None, 'Dioncophyllaceae', '#1030e8'),
    ('DD', ('Drosophyllum', 'Triphyophyllum'), None, '#2050ff', None, None),
    ('Cephalotaceae', ('Cephalotus',), '#f08a00', None, 'Cephalotaceae', '#f08a00'),
    ('Bromeliaceae', ('Brocchinia',), '#b8a030', '#e8c84a', 'Bromeliaceae', '#b8902a'),
    ('Amborellaceae', ('Amborella',), '#000000', None, 'Amborellaceae (outgroup)', '#000000'),
]
BACKBONE = '#6b1500'


def fig3_illustrated(out):
    m = meta()
    tr = Phylo.read('results/kmer/kmer_nj.nwk', 'newick')
    tr.root_with_outgroup('Amborella_trichopoda'); tr.ladderize(reverse=True)
    a_, b_ = tr.root.clades; tot = (a_.branch_length or 0) + (b_.branch_length or 0)
    a_.branch_length = b_.branch_length = tot / 2          # place the root midway on the outgroup branch
    terms = tr.get_terminals(); n = len(terms)
    yy = {t: i for i, t in enumerate(terms)}
    dep = tr.depths(); DMAX = max(dep[t] for t in terms)
    X = {c: dep[c] / DMAX for c in tr.find_clades()}   # phylogram: x proportional to Mash distance
    def sety(c):
        if c in yy: return yy[c]
        v = [sety(ch) for ch in c.clades]; yy[c] = (min(v) + max(v)) / 2; return yy[c]
    sety(tr.root)
    groups = {}
    for key, pre, col, band, flab, fcol in FAM3:
        tips = [t for t in terms if t.name.startswith(pre)]
        groups[key] = (tr.common_ancestor(tips), tips, col, band, flab, fcol)
    colour = {}
    for key, (mrca, tips, col, *_ ) in groups.items():
        if col:
            for c in mrca.find_clades():
                colour[id(c)] = col
    ndh = pd.read_csv('results/ndh_status_curated.csv')
    ndh['taxon'] = ndh.species.str.replace(' x ', '_x_').str.replace(' ', '_'); ndh = ndh.set_index('taxon')

    W_in, H_in = 9.2, 6.6
    fig = plt.figure(figsize=(W_in, H_in))
    ax = fig.add_axes([0.0, 0.06, 1.0, 0.88])
    XMAX = 2.95
    ax.set_xlim(-0.02, XMAX); ax.set_ylim(n + 2.4, -3.6); ax.axis('off')
    ux = (XMAX + 0.02) / W_in; uy = (n + 6.0) / (H_in * 0.88)

    def picture(name, xc, yc, h_rows, z=3):
        img = plt.imread(f'figures/plants/{name}.png')
        h, w = img.shape[:2]
        hd = h_rows; wd = (hd / uy) * (w / h) * ux
        ax.imshow(img, extent=(xc - wd / 2, xc + wd / 2, yc + hd / 2, yc - hd / 2), zorder=z, interpolation='lanczos', aspect='auto')
        return wd

    # clade bands (horizontal gradient)
    for key, (mrca, tips, col, band, *_ ) in groups.items():
        if not band: continue
        y0 = min(yy[t] for t in tips) - 0.5; y1 = max(yy[t] for t in tips) + 0.5
        x0 = X[mrca] - 0.02 if key != 'Lentibulariaceae' else X[mrca] - 0.02
        grad = np.ones((1, 256, 4)); grad[..., :3] = to_rgb(band); grad[0, :, 3] = np.linspace(0, 0.85, 256) ** 1.4
        ax.imshow(grad, extent=(x0, 1.0, y1, y0), aspect='auto', zorder=0, interpolation='bilinear')
    # branches
    def draw(c):
        for ch in c.clades:
            col = colour.get(id(ch), BACKBONE)
            vcol = colour.get(id(c), BACKBONE)
            ax.plot([X[c], X[c]], [yy[c], yy[ch]], color=vcol, lw=0.9, solid_capstyle='butt', zorder=2)
            ax.plot([X[c], X[ch]], [yy[ch], yy[ch]], color=col, lw=0.9, zorder=2)
            draw(ch)
        if not c.is_terminal():
            ax.plot(X[c], yy[c], 'o', ms=1.6, color='#111111', zorder=3)
    draw(tr.root)
    ax.plot([X[tr.root] - 0.02, X[tr.root]], [yy[tr.root]] * 2, color='#111111', lw=0.9)
    sb = 0.05 / DMAX
    ax.plot([0.0, sb], [n + 1.2] * 2, color='#111111', lw=0.9)
    ax.text(sb / 2, n + 1.6, '0.05', ha='center', va='top', fontsize=5.5, color=INK2)
    # tip labels
    for t in terms:
        col = colour.get(id(t), '#000000')
        if X[t] < 0.995:
            ax.plot([X[t], 1.0], [yy[t]] * 2, color='#c9c8c3', lw=0.45, ls=(0, (1, 1.5)), zorder=1)
        ax.add_patch(Rectangle((1.004, yy[t] - 0.33), 0.016, 0.66, color=ORDER_COL[m.loc[t.name, 'order']], lw=0, zorder=3))
        ax.text(1.026, yy[t], pretty(t.name), va='center', fontsize=5.0, color=col, style='italic')
    # ndh matrix
    mx0, cw = 1.43, 0.028
    sym = {'F': dict(fc='#1a1a1a', ec='#1a1a1a'), 'P': dict(fc='white', ec='#1a1a1a'), '-': dict(fc='#ecebe6', ec='#ecebe6')}
    for t in terms:
        r = ndh.loc[t.name]
        for j, g in enumerate(NDH):
            ax.add_patch(Rectangle((mx0 + j * cw + 0.003, yy[t] - 0.36), cw - 0.006, 0.72, lw=0.5, zorder=3, **sym[r[g]]))
    for j, g in enumerate(NDH):
        ax.text(mx0 + j * cw + cw / 2, -0.9, g, rotation=90, ha='center', va='bottom', fontsize=5.0, style='italic')
    ax.text(mx0 + 11 * cw / 2, -3.3, r'Plastid $\it{ndh}$ genes', ha='center', va='bottom', fontsize=6.5, weight='bold')
    # plastome size
    sx0, sw = mx0 + 11 * cw + 0.06, 0.17
    for t in terms:
        L = m.loc[t.name, 'length_bp'] / 1000
        ax.add_patch(Rectangle((sx0, yy[t] - 0.3), sw * (L - 100) / 100, 0.6, fc='#86b6ef', lw=0, zorder=3))
    ax.text(sx0 + sw / 2, -3.3, 'Plastome size', ha='center', va='bottom', fontsize=6.5, weight='bold')
    for v in (100, 150, 200):
        xv = sx0 + sw * (v - 100) / 100
        ax.plot([xv, xv], [n - 0.3, n - 0.1], color=INK2, lw=0.5)
        ax.text(xv, n + 0.4, f'{v}', ha='center', va='top', fontsize=4.6, color=INK2)
    ax.text(sx0 + sw / 2, n + 1.3, 'Plastome size (kb)', ha='center', va='top', fontsize=5.5, color=INK)
    # family brackets and labels
    bx = sx0 + sw + 0.04
    for key, (mrca, tips, col, band, flab, fcol) in groups.items():
        if not flab: continue
        y0 = min(yy[t] for t in tips); y1 = max(yy[t] for t in tips)
        if y1 > y0:
            ax.plot([bx, bx + 0.012, bx + 0.012, bx], [y0 - 0.3, y0 - 0.3, y1 + 0.3, y1 + 0.3], color='#333333', lw=0.6)
        else:
            ax.annotate('', xy=(bx + 0.02, y0), xytext=(bx - 0.015, y0), arrowprops=dict(arrowstyle='-|>', lw=0.6, color='#333333', mutation_scale=5))
        ax.text(bx + 0.03, (y0 + y1) / 2, flab, va='center', fontsize=7.2 if y1 > y0 else 6.2, weight='bold', color=fcol)
    ax.text(XMAX - 0.01, -3.3, 'Major carnivorous plant families', ha='right', va='bottom', fontsize=7, weight='bold')
    # plant drawings: right column next to their families
    px = bx + 0.43
    def yof(prefix):
        ts = [t for t in terms if t.name.startswith(prefix)]
        return (min(yy[t] for t in ts) + max(yy[t] for t in ts)) / 2
    picture('genlisea', px, yof('Genlisea'), 6.5)
    picture('pinguicula_a', px + 0.02, yof('Pinguicula') - 0.3, 4.6)
    picture('pinguicula_b', px + 0.16, yof('Pinguicula') - 0.8, 3.8)
    picture('utricularia', px + 0.05, yof('Utricularia') + 1.0, 9.0)
    picture('sarracenia', px + 0.33, yof('Sarracenia') + 0.3, 5.2)
    picture('roridula', px + 0.14, yof('Roridula') + 0.8, 3.6)
    picture('nepenthes', px + 0.30, yof('Nepenthes'), 4.6)
    picture('drosophyllum', px + 0.09, yof('Drosophyllum') + 0.4, 4.4)
    picture('triphyophyllum', px + 0.20, yof('Triphyophyllum') + 0.6, 4.0)
    picture('brocchinia', px + 0.30, yof('Brocchinia') + 0.2, 4.4)
    picture('cephalotus', px + 0.42, yof('Cephalotus') - 0.6, 3.0)
    picture('amborella', px + 0.12, yof('Amborella') + 0.4, 3.0)
    # drawings within the tree area, beside their clades
    dro = [t for t in terms if t.name.startswith(('Drosera', 'Dionaea', 'Aldrovanda'))]
    yd = (min(yy[t] for t in dro) + max(yy[t] for t in dro)) / 2
    picture('aldrovanda', 0.06, yd - 0.5, 3.8); picture('drosera', 0.135, yd - 0.3, 4.2); picture('dionaea', 0.20, yd - 0.3, 4.2)
    lent = [t for t in terms if t.name.startswith(('Genlisea', 'Pinguicula', 'Utricularia'))]
    yl0 = min(yy[t] for t in lent)
    picture('genlisea_root', 0.30, yl0 + 4.5, 8.0, z=1)
    picture('pinguicula_flower', 0.30, yl0 + 13.0, 8.5, z=1)
    ax.set_xlim(-0.02, XMAX); ax.set_ylim(n + 2.4, -3.6); ax.set_aspect('auto')
    # legend
    hs = [Patch(fc='#1a1a1a', ec='#1a1a1a', label='functional'), Patch(fc='white', ec='#1a1a1a', label='pseudogene'),
          Patch(fc='#ecebe6', ec='#ecebe6', label='absent')] + [Patch(fc=c, label=o) for o, c in ORDER_COL.items()]
    ax.legend(handles=hs, loc='upper center', bbox_to_anchor=(0.5, -0.02), ncol=9, frameon=False, fontsize=6,
              handlelength=1.0, columnspacing=1.2)
    fig.savefig(out, bbox_inches='tight', dpi=300); plt.close(fig)


# ---------------- Figure 5, Nature-style multi-panel codon-usage figure ----------------
from matplotlib.patches import Ellipse
from matplotlib.gridspec import GridSpec
FAMCOL = {'Lentibulariaceae': '#2e9e45', 'Sarraceniaceae': '#4a4ad0', 'Roridulaceae': '#d0137f', 'Droseraceae': '#17a2b8',
          'Nepenthaceae': '#8b1a1a', 'Drosophyllaceae': '#1f3fd6', 'Dioncophyllaceae': '#6f86ff', 'Cephalotaceae': '#f08a00',
          'Bromeliaceae': '#b8962a', 'Amborellaceae': '#555555'}


def conf_ellipse(ax, x, y, col, n_std=2.4477):  # 95% for 2 d.f.
    cov = np.cov(x, y); vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]; vals, vecs = vals[order], vecs[:, order]
    ang = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
    w, h = 2 * n_std * np.sqrt(vals)
    ax.add_patch(Ellipse((x.mean(), y.mean()), w, h, angle=ang, fc=col, ec=col, alpha=0.12, lw=0.8, zorder=1))
    ax.add_patch(Ellipse((x.mean(), y.mean()), w, h, angle=ang, fc='none', ec=col, alpha=0.8, lw=0.7, zorder=1))


def fig5_nature(out):
    rc = {'font.family': 'Liberation Sans', 'font.size': 7, 'axes.linewidth': 0.5, 'xtick.major.width': 0.5,
          'ytick.major.width': 0.5, 'xtick.major.size': 2, 'ytick.major.size': 2}
    with plt.rc_context(rc):
        m = meta()
        r = pd.read_csv('results/codon/rscu.csv', index_col=0)
        r.index = r.index.str.replace(' x ', '_x_').str.replace(' ', '_')
        tr = Phylo.read('results/phylogeny/cds_ml.treefile', 'newick'); tr.root_with_outgroup('Amborella_trichopoda'); tr.ladderize(reverse=True)
        order = [t.name for t in tr.get_terminals() if t.name in r.index]
        B_ = 'TCAG'; AA = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
        code = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate(B_) for j, b in enumerate(B_) for k, c in enumerate(B_)}
        aa3 = dict(F='Phe', L='Leu', S='Ser', Y='Tyr', C='Cys', P='Pro', H='His', Q='Gln', R='Arg', I='Ile', T='Thr', N='Asn',
                   K='Lys', V='Val', A='Ala', D='Asp', E='Glu', G='Gly')
        cods = sorted(r.columns, key=lambda c: (aa3[code[c]], c[2] in 'GC', c))
        M = r.loc[order, cods].T.values            # codons x species
        fig = plt.figure(figsize=(7.2, 7.9))
        gs = GridSpec(2, 3, figure=fig, height_ratios=[1.3, 1], hspace=0.5, wspace=0.5,
                      left=0.085, right=0.975, top=0.965, bottom=0.12)
        # ---- a: heatmap
        gsa = gs[0, :].subgridspec(1, 2, width_ratios=[1, 0.012], wspace=0.01)
        ax = fig.add_subplot(gsa[0, 0])
        cmap = LinearSegmentedColormap.from_list('div', ['#1c5cab', '#86b6ef', '#f4f3ef', '#f2a07b', '#c7381b'])
        vmax = float(np.ceil(np.nanmax(M) * 2) / 2)
        im = ax.imshow(M, aspect='auto', cmap=cmap, norm=TwoSlopeNorm(vmin=0, vcenter=1, vmax=vmax), interpolation='nearest')
        ax.set_yticks(range(len(cods))); ax.set_yticklabels([c.replace('T', 'U') for c in cods], fontsize=4.3, family='DejaVu Sans Mono')
        ax.set_xticks(range(len(order))); ax.set_xticklabels([pretty(t) for t in order], rotation=90, fontsize=4.4, style='italic')
        ax.tick_params(length=0, pad=1)
        for k in range(1, len(cods)):
            if code[cods[k]] != code[cods[k - 1]]:
                ax.axhline(k - 0.5, color='white', lw=0.8)
        # amino-acid labels and third-position bar on the left
        prev = None; start = 0
        for k, c in enumerate(cods + ['END']):
            a = code.get(c)
            if a != prev and prev is not None:
                ax.text(-6.6, (start + k - 1) / 2, aa3[prev], ha='right', va='center', fontsize=5, clip_on=False)
                start = k
            prev = a
        for k, c in enumerate(cods):
            ax.add_patch(Rectangle((-5.4, k - 0.5), 0.9, 1, color='#3a3a3a' if c[2] in 'AT' else '#c9c8c3', clip_on=False, lw=0))
        # family and order bars on top
        for jx, t in enumerate(order):
            ax.add_patch(Rectangle((jx - 0.5, -3.4), 1, 1.3, color=FAMCOL[m.loc[t, 'family']], clip_on=False, lw=0))
            ax.add_patch(Rectangle((jx - 0.5, -1.9), 1, 1.3, color=ORDER_COL[m.loc[t, 'order']], clip_on=False, lw=0))
        ax.text(-0.9, -2.75, 'Family', ha='right', va='center', fontsize=5.5, clip_on=False)
        ax.text(-0.9, -1.25, 'Order', ha='right', va='center', fontsize=5.5, clip_on=False)
        ax.set_ylim(len(cods) - 0.5, -3.6)
        for sp in ax.spines.values(): sp.set_visible(False)
        cax = fig.add_subplot(gsa[0, 1])
        cb = fig.colorbar(im, cax=cax); cb.outline.set_linewidth(0.4); cb.ax.tick_params(labelsize=5.5, length=2)
        cb.set_label('RSCU', fontsize=6)
        # ---- b: PCA with ellipses and loadings
        X = r.loc[order].values; mu = X.mean(0); sd = X.std(0, ddof=1); keep = sd > 0
        Z = (X[:, keep] - mu[keep]) / sd[keep]
        U, S, Vt = np.linalg.svd(Z, full_matrices=False); var = S ** 2 / (S ** 2).sum()
        sc = U * S; load = Vt.T * S / np.sqrt(len(Z) - 1)
        names_k = np.array(r.columns)[keep]
        axb = fig.add_subplot(gs[1, 0])
        fams = pd.Series([m.loc[t, 'family'] for t in order])
        for f in fams.unique():
            idx = np.where(fams == f)[0]
            if len(idx) >= 3:
                conf_ellipse(axb, sc[idx, 0], sc[idx, 1], FAMCOL[f])
            axb.scatter(sc[idx, 0], sc[idx, 1], s=11, color=FAMCOL[f], edgecolor='white', lw=0.4, zorder=3)
        k_ = order.index('Nepenthes_khasiana')
        axb.annotate('N. khasiana', (sc[k_, 0], sc[k_, 1]), xytext=(-22, 16), textcoords='offset points', fontsize=5.5, style='italic',
                     arrowprops=dict(arrowstyle='-', lw=0.4, color=INK2))
        top = np.argsort(-(load[:, 0] ** 2 + load[:, 1] ** 2))[:6]
        f_ = 0.75 * min(np.abs(sc[:, 0]).max() / np.abs(load[top, 0]).max(), np.abs(sc[:, 1]).max() / np.abs(load[top, 1]).max())
        groups_ = []
        for t in top:
            ang = np.degrees(np.arctan2(load[t, 1], load[t, 0]))
            for gq in groups_:
                if abs(((ang - gq['ang']) + 180) % 360 - 180) < 15:
                    gq['m'].append(t); break
            else:
                groups_.append({'ang': ang, 'm': [t]})
        for t in top:
            axb.annotate('', xy=(load[t, 0] * f_, load[t, 1] * f_), xytext=(0, 0),
                         arrowprops=dict(arrowstyle='-|>', lw=0.5, color='#444444', mutation_scale=5), zorder=4)
        for gq in groups_:
            tt = max(gq['m'], key=lambda t: load[t, 0] ** 2 + load[t, 1] ** 2)
            lx_, ly_ = load[tt, 0] * f_, load[tt, 1] * f_
            ha = 'left' if lx_ > 0 else 'right'
            axb.text(lx_ + (0.4 if lx_ > 0 else -0.4), ly_ + (0.5 if ly_ >= 0 else -0.8), ', '.join(names_k[t].replace('T', 'U') for t in gq['m']),
                     fontsize=5, color='#333333', ha=ha, va='center', family='DejaVu Sans Mono', zorder=4,
                     bbox=dict(boxstyle='square,pad=0.1', fc='white', ec='none', alpha=0.8))
        tipx = np.r_[sc[:, 0], load[top, 0] * f_]; tipy = np.r_[sc[:, 1], load[top, 1] * f_]
        ax0, ax1 = axb.get_xlim(); ay0, ay1 = axb.get_ylim()
        axb.set_xlim(min(ax0, tipx.min() - 2.5), max(ax1, tipx.max() + 6.5)); axb.set_ylim(min(ay0, tipy.min() - 1.5), max(ay1, tipy.max() + 1.8))
        axb.axhline(0, color=GRID, lw=0.4, zorder=0); axb.axvline(0, color=GRID, lw=0.4, zorder=0)
        axb.set_xlabel(f'PC1 ({var[0]:.1%})'); axb.set_ylabel(f'PC2 ({var[1]:.1%})')
        for sp in ('top', 'right'): axb.spines[sp].set_visible(False)
        axb.tick_params(labelsize=6)
        # ---- c: ENC vs GC3s
        e = pd.read_csv('results/codon/enc_gc3.csv'); e['family'] = e.species.map(dict(zip(m.species, m.family)))
        axc = fig.add_subplot(gs[1, 1])
        s_ = np.linspace(0.2, 0.36, 200)
        axc.plot(s_, 2 + s_ + 29 / (s_ ** 2 + (1 - s_) ** 2), color='#333333', lw=0.8, zorder=1)
        xe_ = 0.300; ye_ = 2 + xe_ + 29 / (xe_ ** 2 + (1 - xe_) ** 2)
        axc.annotate('Expected under\nmutation only', xy=(xe_, ye_), xytext=(0.2445, 52.9), fontsize=5.2, color='#333333',
                     va='center', ha='left', arrowprops=dict(arrowstyle='-', lw=0.4, color='#777777'))
        for f, g in e.groupby('family'):
            axc.scatter(g.GC3s, g.ENC, s=11, color=FAMCOL[f], edgecolor='white', lw=0.4, zorder=3,
                        marker='D' if f == 'Amborellaceae' else 'o')
        axc.set_xlim(0.24, 0.31); axc.set_ylim(47.5, 53.5)
        axc.set_xlabel('GC3s'); axc.set_ylabel('ENC')
        for sp in ('top', 'right'): axc.spines[sp].set_visible(False)
        axc.tick_params(labelsize=6)
        # ---- d: dumbbell of mean RSCU, A/U- vs G/C-ending codons, by family
        axd = fig.add_subplot(gs[1, 2])
        at = [c for c in r.columns if c[2] in 'AT']; gc = [c for c in r.columns if c[2] in 'GC']
        fl = [f for f in ['Lentibulariaceae', 'Sarraceniaceae', 'Roridulaceae', 'Droseraceae', 'Nepenthaceae', 'Drosophyllaceae',
                          'Dioncophyllaceae', 'Cephalotaceae', 'Bromeliaceae'] if f in set(fams)]
        for i, f in enumerate(fl):
            ts = [t for t in order if m.loc[t, 'family'] == f]
            va = r.loc[ts, at].mean(axis=1).values; vg = r.loc[ts, gc].mean(axis=1).values
            axd.plot([vg.mean(), va.mean()], [i, i], color='#bdbcb6', lw=1.6, zorder=1, solid_capstyle='round')
            axd.scatter(va.mean(), i, s=22, color=FAMCOL[f], edgecolor='white', lw=0.5, zorder=3)
            axd.scatter(vg.mean(), i, s=22, facecolor='white', edgecolor=FAMCOL[f], lw=1.0, zorder=3)
            if len(ts) > 2:
                axd.plot([va.min(), va.max()], [i, i], color=FAMCOL[f], lw=0.5, zorder=2)
                axd.plot([vg.min(), vg.max()], [i, i], color=FAMCOL[f], lw=0.5, zorder=2)
        axd.axvline(1, color='#888888', lw=0.5, ls=(0, (2, 2)))
        axd.set_yticks(range(len(fl))); axd.set_yticklabels(fl, fontsize=5.6); axd.invert_yaxis()
        axd.set_xlim(0.4, 1.6); axd.set_xlabel('Mean RSCU')
        axd.text(0.58, -0.9, 'G/C-ending', ha='center', fontsize=5.4, color='#333333')
        axd.text(1.43, -0.9, 'A/U-ending', ha='center', fontsize=5.4, color='#333333')
        axd.set_ylim(len(fl) - 0.5, -1.3)
        for sp in ('top', 'right'): axd.spines[sp].set_visible(False)
        axd.tick_params(labelsize=6)
        for axx, lab in ((ax, 'A'), (axb, 'B'), (axc, 'C'), (axd, 'D')):
            p = axx.get_position()
            fig.text(0.008 if lab in 'AB' else p.x0 - 0.075, (p.y1 + 0.035) if lab == 'A' else p.y1 + 0.02, lab, fontsize=9, weight='bold')
        # shared family legend
        hs = [Line2D([], [], ls='', marker='o', color=FAMCOL[f], markeredgecolor='white', markersize=4.5, label=f) for f in fl]
        hs += [Line2D([], [], ls='', marker='D', color=FAMCOL['Amborellaceae'], markeredgecolor='white', markersize=4, label='Amborellaceae (outgroup)')]
        hs += [Patch(fc=c, label=o) for o, c in ORDER_COL.items() if o != 'Amborellales']
        hs += [Patch(fc='#3a3a3a', label='A/U-ending codon'), Patch(fc='#c9c8c3', label='G/C-ending codon')]
        fig.legend(handles=hs, loc='lower center', ncol=6, frameon=False, fontsize=5.6, bbox_to_anchor=(0.5, 0.0),
                   handletextpad=0.3, columnspacing=1.0)
        fig.savefig(out, dpi=300, bbox_inches='tight'); plt.close(fig)


def fig6_nature(out):
    """Figure 6: whole-plastome ANI. A, clustered matrix (lower: ANI; upper: aligned fraction);
    B, ANI by taxonomic level; C, ANI vs ML patristic distance; D, aligned fraction vs ANI."""
    from matplotlib.gridspec import GridSpec
    from matplotlib.colors import Normalize
    import itertools
    rc = {'font.family': 'Liberation Sans', 'font.size': 7, 'axes.linewidth': 0.5, 'xtick.major.width': 0.5,
          'ytick.major.width': 0.5, 'xtick.major.size': 2, 'ytick.major.size': 2}
    with plt.rc_context(rc):
        m = meta()
        A = pd.read_csv('results/ani/ani_matrix.tsv', sep='\t', index_col=0)
        F = pd.read_csv('results/ani/alignment_fraction.tsv', sep='\t', index_col=0).loc[A.index, A.index]
        Asym = (A.values + A.values.T) / 2; Fsym = (F.values + F.values.T) / 2
        D = 100 - Asym; np.fill_diagonal(D, 0)
        Z = linkage(squareform(D, checks=False), method='average')
        o = leaves_list(Z); names = list(A.index[o]); n = len(names)
        Ma = Asym[np.ix_(o, o)]; Mf = Fsym[np.ix_(o, o)]
        fam = [m.loc[t, 'family'] for t in names]; ordr = [m.loc[t, 'order'] for t in names]
        # comparison classes
        tr = Phylo.read('results/phylogeny/cds_ml.treefile', 'newick')
        tips = {t.name: t for t in tr.get_terminals()}
        rows = []
        for i, j in itertools.combinations(range(len(A.index)), 2):
            a, b = A.index[i], A.index[j]
            if 'Amborella' in a or 'Amborella' in b: cls = 'Outgroup'
            elif a.split('_')[0] == b.split('_')[0]: cls = 'Within genus'
            elif m.loc[a, 'family'] == m.loc[b, 'family']: cls = 'Within family'
            elif m.loc[a, 'order'] == m.loc[b, 'order']: cls = 'Within order'
            else: cls = 'Between orders'
            rows.append(dict(a=a, b=b, ani=Asym[i, j], af=Fsym[i, j], cls=cls, pd=tr.distance(tips[a], tips[b])))
        P = pd.DataFrame(rows)
        P.to_csv('results/ani/pairwise_ani_classes.csv', index=False)
        CLS = ['Within genus', 'Within family', 'Within order', 'Between orders', 'Outgroup']
        CC = {'Within genus': '#0d366b', 'Within family': '#2a78d6', 'Within order': '#1baf7a',
              'Between orders': '#eda100', 'Outgroup': '#8a8984'}

        fig = plt.figure(figsize=(7.2, 9.3))
        gsT = GridSpec(1, 1, figure=fig, left=0.22, right=0.90, top=0.985, bottom=0.365)
        gs = GridSpec(1, 3, figure=fig, width_ratios=[1.2, 1, 1], wspace=0.42, left=0.08, right=0.985, top=0.235, bottom=0.07)
        # ---- A: matrix
        axA = fig.add_subplot(gsT[0, 0])
        cA = LinearSegmentedColormap.from_list('ani', ['#f7f4ea', '#c9e3d6', '#6fb7c2', '#2a78d6', '#1b3f8c', '#0b1a3d'])
        cF = LinearSegmentedColormap.from_list('af', ['#fbf7f2', '#f6d7bf', '#eda16b', '#c75b26', '#6e2a0c'])
        low = np.where(np.tril(np.ones((n, n), bool), -1), Ma, np.nan)
        up = np.where(np.triu(np.ones((n, n), bool), 1), Mf * 100, np.nan)
        imA = axA.imshow(low, cmap=cA, vmin=80, vmax=100, interpolation='nearest', aspect='equal')
        imF = axA.imshow(up, cmap=cF, vmin=20, vmax=100, interpolation='nearest', aspect='equal')
        for k in range(n):
            axA.add_patch(Rectangle((k - 0.5, k - 0.5), 1, 1, fc='#ffffff', ec='none', zorder=3))
        # family blocks on the diagonal
        k = 0
        while k < n:
            e = k
            while e + 1 < n and fam[e + 1] == fam[k]: e += 1
            axA.add_patch(Rectangle((k - 0.5, k - 0.5), e - k + 1, e - k + 1, fc='none', ec=INK, lw=0.6, zorder=4))
            k = e + 1
        axA.set_xlim(-0.5, n - 0.5); axA.set_ylim(n - 0.5, -0.5)
        axA.set_yticks(range(n)); axA.set_yticklabels([pretty(t) for t in names], fontsize=5.0, style='italic')
        for tl, f in zip(axA.get_yticklabels(), fam): tl.set_color(FAMCOL[f])
        axA.set_xticks([]); axA.tick_params(length=0)
        for sp in axA.spines.values(): sp.set_visible(False)
        pA = axA.get_position()
        # family/order strips beneath matrix
        h = pA.height / n
        axS = fig.add_axes([pA.x0, pA.y0 - 2.6 * h, pA.width, 2.2 * h])
        axS.imshow(np.array([[matplotlib.colors.to_rgb(FAMCOL[f]) for f in fam],
                             [matplotlib.colors.to_rgb(ORDER_COL[x]) for x in ordr]]), aspect='auto', interpolation='nearest')
        axS.set_yticks([0, 1]); axS.set_yticklabels(['Family', 'Order'], fontsize=5.2); axS.set_xticks([])
        axS.tick_params(length=0); [sp.set_visible(False) for sp in axS.spines.values()]
        # dendrogram on the right of the matrix
        axD = fig.add_axes([pA.x1 + 0.004, pA.y0, 0.07, pA.height])
        dendrogram(Z, ax=axD, orientation='right', color_threshold=0, above_threshold_color=INK2, no_labels=True)
        axD.invert_yaxis(); axD.axis('off')
        for ln in axD.collections: ln.set_linewidth(0.5)
        # colour bars
        cx = pA.x0 - 0.205
        c1 = fig.add_axes([pA.x0 + 0.02, pA.y0 - 0.058, 0.22, 0.008]); cb1 = fig.colorbar(imA, cax=c1, orientation='horizontal')
        cb1.set_label('ANI, lower triangle (%)', fontsize=5.6, labelpad=1.5); cb1.ax.tick_params(labelsize=5.2, length=1.5); cb1.outline.set_linewidth(0.4)
        c2 = fig.add_axes([pA.x1 - 0.24, pA.y0 - 0.058, 0.22, 0.008]); cb2 = fig.colorbar(imF, cax=c2, orientation='horizontal')
        cb2.set_label('Aligned fraction, upper triangle (%)', fontsize=5.6, labelpad=1.5); cb2.ax.tick_params(labelsize=5.2, length=1.5); cb2.outline.set_linewidth(0.4)

        rng = np.random.default_rng(1)
        # ---- B: ANI by level
        axB = fig.add_subplot(gs[0, 0])
        for i_, c in enumerate(CLS):
            v = P.loc[P.cls == c, 'ani'].values
            axB.scatter(i_ + rng.uniform(-0.28, 0.28, len(v)), v, s=2.2, color=CC[c], alpha=0.45, lw=0, zorder=2, rasterized=True)
            q1, q2, q3 = np.percentile(v, [25, 50, 75])
            axB.add_patch(Rectangle((i_ - 0.32, q1), 0.64, q3 - q1, fc='none', ec=INK, lw=0.6, zorder=3))
            axB.plot([i_ - 0.32, i_ + 0.32], [q2, q2], color=INK, lw=1.0, zorder=3)
            axB.text(i_, 101.0, f'{len(v)}', ha='center', fontsize=4.8, color=INK2)
        axB.text(-0.95, 101.0, 'n =', ha='right', fontsize=4.8, color=INK2)
        axB.set_xticks(range(len(CLS))); axB.set_xticklabels(['Same\ngenus', 'Same\nfamily', 'Same\norder', 'Different\norders', 'vs\noutgroup'], fontsize=5.2)
        axB.set_ylabel('ANI (%)'); axB.set_ylim(78, 102)
        axB.set_xlim(-0.6, len(CLS) - 0.4)
        # ---- C: ANI vs patristic distance
        axC = fig.add_subplot(gs[0, 1])
        for c in CLS:
            q = P[P.cls == c]
            axC.scatter(q.pd, q.ani, s=3, color=CC[c], alpha=0.6, lw=0, label=c, rasterized=True)
        xs = np.sort(P.pd.values); bins = np.quantile(xs, np.linspace(0, 1, 13))
        mid = [(P.pd[(P.pd >= a) & (P.pd <= b)]).median() for a, b in zip(bins[:-1], bins[1:])]
        med = [(P.ani[(P.pd >= a) & (P.pd <= b)]).median() for a, b in zip(bins[:-1], bins[1:])]
        axC.plot(mid, med, color=INK, lw=0.9, zorder=4)
        rho = pd.Series(P.pd).corr(pd.Series(P.ani), method='spearman')
        axC.text(0.97, 0.96, f'Spearman ρ = {rho:.2f}', transform=axC.transAxes, ha='right', va='top', fontsize=5.4)
        axC.set_xlabel('Patristic distance (subst./site)'); axC.set_ylabel('ANI (%)'); axC.set_ylim(78, 102)
        # ---- D: aligned fraction vs ANI
        axDd = fig.add_subplot(gs[0, 2])
        for c in CLS:
            q = P[P.cls == c]
            axDd.scatter(q.ani, q.af * 100, s=3, color=CC[c], alpha=0.6, lw=0, rasterized=True)
        rho2 = pd.Series(P.ani).corr(pd.Series(P.af), method='spearman')
        axDd.text(0.97, 0.04, f'Spearman ρ = {rho2:.2f}', transform=axDd.transAxes, ha='right', va='bottom', fontsize=5.4)
        axDd.set_xlabel('ANI (%)'); axDd.set_ylabel('Aligned fraction (%)'); axDd.set_xlim(78, 101); axDd.set_ylim(0, 102)
        for axx in (axB, axC, axDd):
            for sp in ('top', 'right'): axx.spines[sp].set_visible(False)
            axx.tick_params(labelsize=6)
        LBL = dict(zip(CLS, ['Same genus', 'Same family', 'Same order', 'Different orders', 'Versus outgroup']))
        hs = [Line2D([], [], ls='', marker='o', color=CC[c], markersize=4, label=LBL[c]) for c in CLS]
        fig.legend(handles=hs, loc='lower center', ncol=5, frameon=False, fontsize=5.8, bbox_to_anchor=(0.53, 0.0),
                   handletextpad=0.2, columnspacing=1.2)
        hf = [Patch(fc=FAMCOL[f], label=f) for f in dict.fromkeys(fam)] + [Patch(fc=ORDER_COL[x], label=x) for x in dict.fromkeys(ordr)]
        fig.legend(handles=hf, loc='upper center', ncol=8, frameon=False, fontsize=5.4, bbox_to_anchor=(0.53, pA.y0 - 0.085),
                   handlelength=1.0, handleheight=0.8, handletextpad=0.4, labelspacing=0.3, columnspacing=1.0)
        fig.text(0.015, pA.y1 + 0.004, 'A', fontsize=9, weight='bold')
        for axx, lab in ((axB, 'B'), (axC, 'C'), (axDd, 'D')):
            p_ = axx.get_position(); fig.text(p_.x0 - 0.065, p_.y1 + 0.022, lab, fontsize=9, weight='bold')
        fig.savefig(out, dpi=300, bbox_inches='tight'); plt.close(fig)
        print(P.groupby('cls').ani.describe().round(2)); print('rho_pd', round(rho, 3), 'rho_af', round(rho2, 3))
        print(P.groupby('cls').af.median().round(3))



if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--outdir', default='figures'); ap.add_argument('--only', default='')
    a = ap.parse_args(); os.makedirs(a.outdir, exist_ok=True)
    todo = a.only.split(',') if a.only else ['3', '5', '6', 'S1', '2A']
    if '3' in todo: fig3_illustrated(os.path.join(a.outdir, 'Figure3.png'))
    if '3plain' in todo: fig3(os.path.join(a.outdir, 'Figure3_plain.png'))
    if '5' in todo: fig5_nature(os.path.join(a.outdir, 'Figure5.png'))
    if '5plain' in todo: fig5(os.path.join(a.outdir, 'Figure5_plain.png'))
    if '6' in todo: fig6_nature(os.path.join(a.outdir, 'Figure6.png'))
    if '6plain' in todo: fig6(os.path.join(a.outdir, 'Figure6_plain.png'))
    if 'S1' in todo: figS1(os.path.join(a.outdir, 'FigureS1.png'))
    if '2A' in todo: fig2a_classic(os.path.join(a.outdir, 'Figure2A.png'))
    if '2Aflat' in todo: fig2a(os.path.join(a.outdir, 'Figure2A_alltips.png'))


