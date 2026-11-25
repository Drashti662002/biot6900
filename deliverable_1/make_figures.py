"""
make_figures.py

Re-plots existing Lab 2 and Lab 4 results for Deliverable 1. No new analysis:
every value comes from an existing lab output (see provenance.md).

    Figure 1  A: RNA vs protein log2FC for the 4,851 joined genes (Lab 2)
              B: three-set Venn of FDR < 0.05 support per layer (Lab 2 q-values)
    Figure 2  AlphaFold models of CA2 and CR2 coloured by pLDDT, with Lab 4's
              fpocket pockets marked, plus per-residue pLDDT traces

Run from the repo root (needs PyMOL.app for the structure renders):

    python deliverable_1/make_figures.py
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "deliverable_1/figures"
STRUCT = ROOT / "shared/module4_structures"
PYMOL = "/Applications/PyMOL.app/Contents/MacOS/PyMOL"
OUT.mkdir(parents=True, exist_ok=True)

# Okabe-Ito colours; CVD separation checked (all pairs dE > 8 under protan/deutan/tritan)
BLUE, VERMILLION, ORANGE, GREEN, PURPLE = "#0072B2", "#D55E00", "#E69F00", "#009E73", "#CC79A7"
INK, MUTED, GRID = "#222222", "#666666", "#DDDDDD"
# AlphaFold pLDDT bands (standard AFDB colours)
AF_BANDS = [(90, 100, "#0053D6", "Very high (>90)"), (70, 90, "#65CBF3", "Confident (70-90)"),
            (50, 70, "#FFDB13", "Low (50-70)"), (0, 50, "#FF7D45", "Very low (<50)")]

plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8, "xtick.labelsize": 7.5,
                     "ytick.labelsize": 7.5, "legend.fontsize": 7.5, "font.family": "DejaVu Sans",
                     "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "pdf.fonttype": 42})


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- Figure 1

def figure1():
    t = pd.read_csv(ROOT / "modules/module2/targets_coad.csv")
    lab3 = pd.read_csv(ROOT / "modules/module3/targets_coad_w3.csv")["gene"].tolist()
    both_sig = (t.rna_q < 0.05) & (t.prot_q < 0.05)

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(6.6, 3.3), gridspec_kw={"width_ratios": [1.25, 1]})

    # --- A: scatter
    ns = t[~both_sig]
    ax.scatter(ns.rna_lfc, ns.prot_lfc, s=3, c="#BBBBBB", lw=0, label=f"Not significant on both (n = {len(ns):,})",
               rasterized=True)
    for conc, col, mk, lab in [(True, BLUE, "o", "Concordant"), (False, VERMILLION, "^", "Discordant")]:
        d = t[both_sig & (t.concordant == conc)]
        ax.scatter(d.rna_lfc, d.prot_lfc, s=4 if conc else 5, c=col, marker=mk, lw=0, alpha=0.75,
                   label=f"{lab}, both FDR < 0.05 (n = {len(d):,})", rasterized=True)
    top = t[t.gene.isin(lab3)]
    ax.scatter(top.rna_lfc, top.prot_lfc, s=22, facecolor="none", edgecolor=INK, lw=0.8,
               label="Top 15 carried to triage")
    for g, dx, dy in [("CA2", 0.25, -0.05), ("CR2", 0.3, 0.25), ("DES", 0.45, -0.25), ("SOX9", 0.15, 0.25),
                      ("COL5A1", 0.15, -0.35)]:
        r = t.loc[t.gene == g].iloc[0]
        ax.annotate(g, (r.rna_lfc, r.prot_lfc), xytext=(r.rna_lfc + dx, r.prot_lfc + dy), fontsize=7.5,
                    color=INK, fontweight="bold" if g in ("CA2", "CR2") else "normal",
                    arrowprops=dict(arrowstyle="-", lw=0.5, color=MUTED))
    ax.axhline(0, color=MUTED, lw=0.6); ax.axvline(0, color=MUTED, lw=0.6)
    ax.set_xlabel("RNA log2 fold change, tumour vs normal (TCGA-COAD)")
    ax.set_ylabel("Protein log2 fold change, tumour vs normal (CPTAC)")
    ax.set_title("A  RNA-protein concordance, 4,851 genes", loc="left", fontweight="bold")
    ax.legend(loc="upper left", frameon=False, fontsize=6.5, handletextpad=0.3, borderaxespad=0.1, markerscale=1.6)
    ax.set_xlim(-7.6, 6.2); ax.set_ylim(-3.2, 4.4)
    for s in ("top", "right"): ax.spines[s].set_visible(False)

    # --- B: Venn (FDR < 0.05 per layer)
    R, P, G = t.rna_q < 0.05, t.prot_q < 0.05, t.gen_q < 0.05
    n = {"r": (R & ~P & ~G).sum(), "p": (~R & P & ~G).sum(), "g": (~R & ~P & G).sum(),
         "rp": (R & P & ~G).sum(), "rg": (R & ~P & G).sum(), "pg": (~R & P & G).sum(), "rpg": (R & P & G).sum(),
         "none": (~R & ~P & ~G).sum()}
    centers = {"r": (-0.55, 0.32), "p": (0.55, 0.32), "g": (0.0, -0.58)}
    for k, col, lab, cnt in [("r", BLUE, "Transcriptomic", R.sum()), ("p", ORANGE, "Proteomic", P.sum()),
                             ("g", GREEN, "Genomic", G.sum())]:
        bx.add_patch(Circle(centers[k], 1.0, facecolor=col, alpha=0.18, edgecolor=col, lw=1.2))
    bx.text(-1.5, 1.2, f"Transcriptomic\nn = {R.sum():,}", ha="center", color=INK, fontsize=7.5)
    bx.text(1.5, 1.2, f"Proteomic\nn = {P.sum():,}", ha="center", color=INK, fontsize=7.5)
    bx.text(0.0, -1.83, f"Genomic (mutation) n = {G.sum()}", ha="center", color=INK, fontsize=7.5)
    pos = {"r": (-1.05, 0.55), "p": (1.05, 0.55), "g": (0, -1.15), "rp": (0, 0.6), "rg": (-0.62, -0.42),
           "pg": (0.62, -0.42), "rpg": (0, -0.05)}
    for k, (x, y) in pos.items():
        bx.text(x, y, f"{n[k]:,}", ha="center", va="center", fontsize=8.5, color=INK,
                fontweight="bold" if k in ("rp", "rpg") else "normal")
    bx.annotate("All 15 triage\ncandidates", xy=(0.0, 0.9), xytext=(0.0, 1.75), ha="center", fontsize=7,
                color=INK, arrowprops=dict(arrowstyle="->", lw=0.6, color=INK))
    bx.annotate("16 genes, incl. SOX9,\nKRAS, SMAD4;\nnone carried forward", xy=(0.08, -0.12),
                xytext=(1.25, -1.2), fontsize=6.5, color=INK, ha="left",
                arrowprops=dict(arrowstyle="->", lw=0.6, color=INK))
    bx.text(-1.75, -1.45, f"None: {n['none']}", fontsize=7, color=MUTED)
    bx.set_xlim(-2.1, 2.3); bx.set_ylim(-2.0, 2.3); bx.set_aspect("equal"); bx.axis("off")
    bx.set_title("B  Genes significant per layer (FDR < 0.05)", loc="left", fontweight="bold")
    fig.tight_layout(w_pad=0.8)
    save(fig, "figure1_integration_venn")
    return n


# ---------------------------------------------------------------- Figure 2

def ca_plddt(pdb):
    return pd.Series({int(l[22:26]): float(l[60:66]) for l in open(pdb)
                      if l.startswith("ATOM") and l[12:16].strip() == "CA"}).sort_index()


def pocket_res(acc, p):
    f = STRUCT / f"AF-{acc}_out/pockets/pocket{p}_atm.pdb"
    return sorted({int(l[22:26]) for l in open(f) if l.startswith("ATOM")})


def render(acc, pockets, png, ligand=None, turn="", buffer=2, spheres=False):
    """pockets: list of (residues, colour). ligand: (pdb, selection) already in the same frame."""
    model = STRUCT / (f"AF-{acc}_aligned_to_3HS4.pdb" if ligand else f"AF-{acc}.pdb")
    lines = [f"load {model}, m", "hide everything", "show cartoon, m",
             "set cartoon_transparency, 0.0", "bg_color white", "set ray_opaque_background, 0",
             "set ray_shadows, 0", "set antialias, 2", "set cartoon_fancy_helices, 1"]
    for i, (lo, hi, col, _) in enumerate(AF_BANDS):
        lines += [f"set_color band{i}, [{int(col[1:3], 16)}, {int(col[3:5], 16)}, {int(col[5:7], 16)}]",
                  f"color band{i}, m and b > {lo - 0.001} and b < {hi + 0.001}"]
    for k, (res, col) in enumerate(pockets):
        sel = "+".join(map(str, res))
        lines += [f"set_color pk{k}, [{int(col[1:3], 16)}, {int(col[3:5], 16)}, {int(col[5:7], 16)}]",
                  f"select p{k}, m and resi {sel} and not name N+C+O",
                  f"show {'spheres' if spheres else 'sticks'}, p{k}", f"color pk{k}, p{k}",
                  "set stick_radius, 0.28", "set sphere_scale, 0.9"]
    if ligand:
        lines += [f"load {ligand[0]}, xtal", f"create lig, xtal and {ligand[1]}", "delete xtal",
                  "show sticks, lig", "color black, lig and elem C", "set stick_radius, 0.3, lig"]
    lines += ["orient m", turn, f"zoom m, {buffer}", f"png {png}, width=1400, height=1200, dpi=300, ray=1", "quit"]
    script = OUT / "_render.pml"
    script.write_text("\n".join(l for l in lines if l))
    subprocess.run([PYMOL, "-cq", str(script)], check=True, capture_output=True)
    script.unlink()


def figure2():
    lab4 = pd.read_csv(ROOT / "modules/module4/structure_assessment_coad.csv").set_index("gene")
    ca2_p2, ca2_p6 = pocket_res("P00918", 2), pocket_res("P00918", 6)
    cr2_p50 = pocket_res("P20023", 50)
    render("P00918", [(ca2_p2, PURPLE), (ca2_p6, VERMILLION)], OUT / "_ca2.png",
           ligand=(STRUCT / "3HS4.pdb", "resn AZM and resi 701 and chain A"))
    render("P20023", [(cr2_p50, PURPLE)], OUT / "_cr2.png", turn="turn z, 90", buffer=12, spheres=True)

    pl = {"CA2": ca_plddt(STRUCT / "AF-P00918.pdb"), "CR2": ca_plddt(STRUCT / "AF-P20023.pdb")}
    fig = plt.figure(figsize=(6.6, 5.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[2.5, 1], hspace=0.32, wspace=0.12)
    panels = [
        ("CA2", "A  CA2 (P00918)", OUT / "_ca2.png",
         f"Mean pLDDT {lab4.loc['CA2', 'mean_pLDDT']}; one domain\nPocket 2 (drug site): score 0.298, pLDDT 97.9\n"
         "Pocket 6 (top score): score 0.317, pLDDT 98.4",
         [(ca2_p2, PURPLE, "Pocket 2 (drug site)"), (ca2_p6, VERMILLION, "Pocket 6 (top score)")]),
        ("CR2", "B  CR2 (P20023)", OUT / "_cr2.png",
         f"Mean pLDDT {lab4.loc['CR2', 'mean_pLDDT']}; 8 domains, PAE 26.3 Å\n"
         "Pocket 50 (top score): score 0.898, pLDDT 65.2\nBest kept pocket 27: score 0.005",
         [(cr2_p50, PURPLE, "Pocket 50 (top score, rejected)")]),
    ]
    for j, (g, title, img, note, pks) in enumerate(panels):
        a = fig.add_subplot(gs[0, j])
        a.imshow(plt.imread(img)); a.axis("off")
        a.set_title(title, loc="left", fontweight="bold")
        a.text(0.0, -0.02, note, transform=a.transAxes, ha="left", va="top", fontsize=6.8, color=INK)
        b = fig.add_subplot(gs[1, j])
        for lo, hi, col, _ in AF_BANDS:
            b.axhspan(lo, hi, color=col, alpha=0.18, lw=0)
        b.plot(pl[g].index, pl[g].values, color=INK, lw=0.7)
        for res, col, lab in pks:
            b.scatter(res, pl[g].loc[res].values, s=10, color=col, zorder=3, label=lab, edgecolor="white", lw=0.3)
        b.set_ylim(0, 100); b.set_xlim(1, pl[g].index.max())
        b.set_xlabel("Residue (UniProt numbering)"); b.set_ylabel("pLDDT" if j == 0 else "")
        b.axhline(70, color=MUTED, lw=0.6, ls="--")
        b.legend(loc="lower left", frameon=False, fontsize=6.5, handletextpad=0.2, borderaxespad=0.1)
        for s in ("top", "right"): b.spines[s].set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, _, c, _ in AF_BANDS]
    fig.legend(handles, [l for *_, l in AF_BANDS], loc="lower center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, -0.03), fontsize=7, title="AlphaFold pLDDT (cartoon colour)",
               title_fontsize=7)
    save(fig, "figure2_structures")
    for f in ("_ca2.png", "_cr2.png"):
        (OUT / f).unlink()


if __name__ == "__main__":
    counts = figure1()
    print("Venn region counts:", counts)
    figure2()
    print("wrote", sorted(p.name for p in OUT.iterdir()))
