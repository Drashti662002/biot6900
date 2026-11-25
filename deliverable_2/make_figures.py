"""Figures for Portfolio Deliverable 2. Reads existing lab outputs only.

No new analysis. The one recomputation is the Lab 6 PRS (PRS = dosages @ weights, with the
allele flips the notebook applied), which is deterministic and is verified in-script against the
AUCs saved in notebook cell 16. If the check fails the script raises instead of plotting.

Run:  python make_figures.py        (writes figures/*.png at 300 dpi and figures/*.pdf)
"""
import os, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, FancyArrowPatch
import networkx as nx
from sklearn.metrics import roc_auc_score, roc_curve

ROOT = "/Users/drashtibhanushali/biot6900"
M5 = f"{ROOT}/modules/module5/data"
M6 = f"{ROOT}/modules/module 6/part3_crc/data"
M7 = f"{ROOT}/modules/module 7"
OUT = f"{ROOT}/deliverable_2/figures"
os.makedirs(OUT, exist_ok=True)

# Okabe-Ito, colourblind-safe
CB = dict(blue="#0072B2", orange="#E69F00", green="#009E73", red="#D55E00",
          purple="#CC79A7", sky="#56B4E9", yellow="#F0E442", grey="#999999", black="#000000")
plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 300, "savefig.dpi": 300, "pdf.fonttype": 42,
})

def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}/{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {name}.png / .pdf")

def panel(ax, letter):
    ax.text(-0.14, 1.06, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")

# ---------------------------------------------------------------- Figure 1
def figure1():
    d = pd.read_csv(f"{M5}/DPYD.tsv", sep="\t")
    sig = d.clinvar_sig.fillna("").str.lower()
    plp = sig.str.contains("pathogenic") & ~sig.str.contains("conflict") & ~sig.str.contains("benign")
    pos = d[plp].copy()
    SHORT = {"stop_gained": "stop gained", "missense_variant": "missense",
             "splice_donor_variant": "splice donor", "splice_acceptor_variant": "splice acceptor",
             "frameshift_variant": "frameshift"}
    pos["cls"] = pos.consequence.map(SHORT).fillna("no consequence\nannotation")
    order = ["stop gained", "missense", "splice donor", "splice acceptor", "frameshift",
             "no consequence\nannotation"]

    fig = plt.figure(figsize=(7.4, 6.0))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.12, 0.88, 0.74], hspace=0.95)
    axA = fig.add_subplot(gs[0]); axB = fig.add_subplot(gs[1]); axC = fig.add_subplot(gs[2])

    # (A) CADD by consequence class, with an explicit unscored row
    UN = -4.5
    rng = np.random.default_rng(0)
    for i, cls in enumerate(order):
        sub = pos[pos.cls == cls]
        sc, un = sub[sub.cadd_phred.notna()], sub[sub.cadd_phred.isna()]
        if len(sc):
            axA.scatter(i + rng.uniform(-.16, .16, len(sc)), sc.cadd_phred, s=20,
                        color=CB["blue"], alpha=.85, edgecolor="white", linewidth=.4, zorder=3)
        if len(un):
            axA.scatter(i + rng.uniform(-.16, .16, len(un)), np.full(len(un), UN), s=22,
                        color=CB["red"], alpha=.9, marker="X", edgecolor="white",
                        linewidth=.4, zorder=3)
        axA.annotate(f"{len(sc)}/{len(sub)}", (i, 47.5), ha="center", fontsize=7,
                     fontweight="bold", color=CB["blue"] if len(sc) else CB["red"])
    axA.axhspan(-7.2, -2.0, color=CB["red"], alpha=.07, zorder=0)
    for b in (10, 20, 25, 30):
        axA.axhline(b, color=CB["grey"], lw=.5, ls="--", alpha=.55)
    axA.set_xticks(range(len(order))); axA.set_xticklabels(order, fontsize=7)
    axA.set_yticks([UN, 0, 10, 20, 30, 40]); axA.set_yticklabels(["no score", "0", "10", "20", "30", "40"])
    axA.set_ylim(-8, 53); axA.set_xlim(-.6, len(order) - .4)
    axA.set_ylabel("CADD PHRED score")
    axA.set_title("A   CADD coverage of the 105 pathogenic $\\it{DPYD}$ variants "
                  "(fraction scored above each class)", loc="left", fontweight="bold", pad=14)
    axA.legend(handles=[Line2D([], [], marker="o", ls="", ms=4, color=CB["blue"], label="scored by CADD (39)"),
                        Line2D([], [], marker="X", ls="", ms=5, color=CB["red"], label="no CADD score (66)")],
               loc="upper center", frameon=False, ncol=2, bbox_to_anchor=(.5, -0.30), handletextpad=.4)

    # (B) pooled LR+ by band
    bands = ["0-10", "10-20", "20-25", "25-30", "30-50"]
    n, plp_n, ben_n = [76, 34, 17, 6, 34], [2, 1, 8, 4, 34], [74, 33, 9, 2, 0]
    lr = [0.07, 0.07, 2.14, 4.82, np.nan]
    xs = np.arange(len(bands))
    axB.bar(xs[:4], lr[:4], color=[CB["blue"], CB["blue"], CB["orange"], CB["orange"]],
            width=.6, zorder=3)
    axB.bar(xs[4], 8.6, color="none", edgecolor=CB["red"], hatch="///", lw=1.2, width=.6, zorder=3)
    axB.axhline(1, color=CB["black"], lw=.85)
    axB.annotate("LR+ = 1, no information", (-.42, 1.3), fontsize=6.6, color=CB["black"])
    for i, v in enumerate(lr):
        lab = f"{v:.2f}" if not np.isnan(v) else "infinite"
        axB.annotate(lab, (i, (v if i < 4 else 8.6) + .3), ha="center", fontsize=7.5, fontweight="bold")
    axB.annotate("34 pathogenic, 0 benign: an empty cell,\nin the band holding the positives\nthat need no model",
                 (2.55, 9.9), ha="center", va="top", fontsize=6.6, color=CB["red"], fontweight="bold")
    axB.set_xticks(xs)
    axB.set_xticklabels([f"{b}\nn={nn}\n{p} path / {bn} benign" for b, nn, p, bn in zip(bands, n, plp_n, ben_n)],
                        fontsize=6.8)
    axB.set_xlabel("CADD PHRED band", labelpad=1)
    axB.set_ylabel("Positive likelihood\nratio (LR+)")
    axB.set_ylim(0, 10.6); axB.set_xlim(-.6, 4.6)
    axB.set_title("B   Pooled five-gene calibration, 167 labelled variants",
                  loc="left", fontweight="bold", pad=8)

    # (C) tag vs causal, laid out as a labelled strip (no overlapping text)
    rows = [("$\\it{DPYD}$ HapB3 c.1236G>A (rs56038477), synonymous tag", 14.03, "Benign, 2 star",
             CB["blue"], "counted as a benign variant"),
            ("$\\it{DPYD}$ HapB3 c.1129-5923C>G (rs75017182), intronic", 0.378, "drug response, 3 star",
             CB["red"], "the allele that skips exon 11"),
            ("$\\it{CYP2C19}$ *2 (rs4244285), synonymous", 5.686, "drug response, 0 star",
             CB["red"], "creates a cryptic splice donor"),
            ("$\\it{CYP2C19}$ *17 (rs12248560)", 3.031, "Benign, 1 star", CB["blue"], "")]
    axC.axvspan(0, 10, color=CB["blue"], alpha=.07); axC.axvspan(10, 20, color=CB["blue"], alpha=.07)
    for b in (10, 20):
        axC.axvline(b, color=CB["grey"], lw=.5, ls="--", alpha=.6)
    for i, (lab, cadd, cv, col, note) in enumerate(rows):
        y = len(rows) - 1 - i
        axC.scatter([cadd], [y], s=48, color=col, zorder=4, edgecolor="white", linewidth=.5)
        axC.annotate(f"CADD {cadd}", (cadd, y), xytext=(-7, 0), textcoords="offset points",
                     ha="right", va="center", fontsize=6.5, fontweight="bold")
        axC.annotate(f"{lab}  |  ClinVar: {cv}" + (f"  |  {note}" if note else ""),
                     (21.5, y), fontsize=6.6, va="center")
    axC.add_patch(FancyArrowPatch((14.03, 2.78), (0.378, 2.22), arrowstyle="<|-|>",
                                  color=CB["red"], lw=1.0, ls=":", mutation_scale=7,
                                  connectionstyle="arc3,rad=-0.45"))
    axC.annotate("in LD: the tag is scored,\nthe causal allele is not",
                 (7.4, 3.42), ha="center", fontsize=6.5, color=CB["red"], fontweight="bold")
    axC.annotate("bands 0-10 and 10-20 carry LR+ 0.07, which reads as evidence toward benign",
                 (-2, -1.05), ha="left", fontsize=6.4, color=CB["grey"], annotation_clip=False)
    axC.set_yticks([]); axC.set_xlim(-2, 78); axC.set_ylim(-1.45, 4.05)
    axC.set_xticks([0, 5, 10, 15, 20]); axC.set_xlabel("CADD PHRED score", labelpad=1)
    axC.spines["left"].set_visible(False)
    axC.set_title("C   The variants the bands misplace", loc="left", fontweight="bold", pad=18)

    save(fig, "figure1_variant_interpretation")

# ---------------------------------------------------------------- Figure 2
POPS = ["EUR", "AFR", "EAS", "SAS", "AMR"]
PC = {"EUR": CB["blue"], "AFR": CB["red"], "EAS": CB["green"], "SAS": CB["orange"], "AMR": CB["purple"]}
SAVED_AUC = {"EUR": 0.685, "AFR": 0.567, "EAS": 0.689, "SAS": 0.606, "AMR": 0.670}
SAVED_OVERALL = 0.647

def lab6_prs():
    """Reproduce the notebook's PRS exactly (cells 10 and 12) and verify against saved AUCs."""
    meta = pd.read_csv(f"{M6}/meta.csv")
    geno = pd.read_csv(f"{M6}/genotypes.csv")
    score = pd.read_csv(f"{M6}/score.csv")
    variants = pd.read_csv(f"{M6}/variants.csv")
    shared = [v for v in score["variant"] if v in geno.columns]
    G = geno[shared].to_numpy(dtype=float).copy()
    s = score.set_index("variant").loc[shared]; v = variants.set_index("variant").loc[shared]
    eff, ref, alt, w = s["effect_allele"].to_numpy(), v["ref"].to_numpy(), v["alt"].to_numpy(), s["weight"].to_numpy()
    keep, flip = (eff == alt), (eff == ref); drop = ~(keep | flip)
    G[:, flip] = 2 - G[:, flip]
    G, w = G[:, ~drop], w[~drop]
    prs = G @ w
    meta["prs"] = (prs - prs.mean()) / prs.std()
    assert int(flip.sum()) == 28 and int(drop.sum()) == 1 and int((keep | flip).sum()) == 61, \
        f"harmonisation mismatch: used {int((keep|flip).sum())} flipped {int(flip.sum())} dropped {int(drop.sum())}"
    y, z = meta["phenotype"].to_numpy(), meta["prs"].to_numpy()
    got = {p: roc_auc_score(y[(meta.ancestry == p).to_numpy()], z[(meta.ancestry == p).to_numpy()]) for p in POPS}
    for p in POPS:
        assert abs(round(got[p], 3) - SAVED_AUC[p]) < 1e-9, f"{p} AUC {got[p]:.3f} != saved {SAVED_AUC[p]}"
    assert abs(round(roc_auc_score(y, z), 3) - SAVED_OVERALL) < 1e-9
    print("  PRS recomputation verified against saved AUCs (cell 16): all 5 ancestries + overall")
    return meta

# calibration quintiles, notebook cell 19
QUINT = {
    "EUR": ([.13, .22, .29, .36, .51], [.13, .23, .26, .36, .52]),
    "AFR": ([.11, .18, .25, .32, .47], [.16, .26, .31, .26, .34]),
    "EAS": ([.11, .17, .23, .31, .45], [.10, .24, .23, .34, .50]),
    "SAS": ([.15, .23, .30, .38, .53], [.32, .32, .30, .39, .58]),
    "AMR": ([.10, .16, .22, .29, .45], [.11, .20, .25, .29, .46]),
}
FLAG = {"EUR": 12.2, "AFR": 8.4, "EAS": 7.2, "SAS": 14.9, "AMR": 6.3}
OE = {"EUR": 1.00, "AFR": 0.99, "EAS": 1.10, "SAS": 1.20, "AMR": 1.08}
R2 = {"EUR": 0.83, "AFR": 0.45, "EAS": 0.76, "SAS": 0.77, "AMR": 0.74}

def figure2():
    meta = lab6_prs()
    y, z = meta["phenotype"].to_numpy(), meta["prs"].to_numpy()
    fig, ax = plt.subplots(1, 3, figsize=(7.4, 2.95))
    fig.subplots_adjust(wspace=0.52, top=0.80, bottom=0.21)

    # (A) ROC
    ax[0].plot([0, 1], [0, 1], color=CB["grey"], lw=.7, ls=":")
    for p_ in POPS:
        m = (meta.ancestry == p_).to_numpy()
        fpr, tpr, _ = roc_curve(y[m], z[m])
        ax[0].plot(fpr, tpr, color=PC[p_], lw=1.4, label=f"{p_}  {SAVED_AUC[p_]:.3f}")
    ax[0].set_xlabel("False positive rate"); ax[0].set_ylabel("True positive rate")
    ax[0].set_title("A   Discrimination", loc="left", fontweight="bold", pad=6)
    leg = ax[0].legend(frameon=False, loc="lower right", handlelength=.9, fontsize=6.4,
                       labelspacing=.25, borderpad=.1, title="AUC", title_fontsize=6.6)
    leg._legend_box.align = "left"
    ax[0].text(.03, .955, f"overall 0.647\nEUR is the training\nancestry", fontsize=6.3,
               transform=ax[0].transAxes, va="top")
    ax[0].set_xlim(0, 1); ax[0].set_ylim(0, 1)
    ax[0].set_xticks([0, .25, .5, .75, 1]); ax[0].set_yticks([0, .25, .5, .75, 1])

    # (B) calibration
    ax[1].plot([.05, .6], [.05, .6], color=CB["grey"], lw=.7, ls=":")
    for p_ in POPS:
        pr, ob = QUINT[p_]
        ax[1].plot(pr, ob, "o-", color=PC[p_], lw=1.3, ms=3.1, label=f"{p_}  {OE[p_]:.2f}")
    pr, ob = QUINT["AFR"]
    ax[1].annotate("AFR: flat.\npredicted .11 to .47,\nobserved .16 to .34", (pr[4], ob[4]),
                   xytext=(.255, .50), fontsize=6.2, color=CB["red"], fontweight="bold",
                   ha="center", arrowprops=dict(arrowstyle="->", color=CB["red"], lw=.75))
    pr, ob = QUINT["SAS"]
    ax[1].annotate("SAS Q1: observed .32\nvs predicted .15 (2.1x)", (pr[0], ob[0]),
                   xytext=(.145, .615), fontsize=6.2, color=CB["orange"], fontweight="bold",
                   ha="left", arrowprops=dict(arrowstyle="->", color=CB["orange"], lw=.75))
    ax[1].set_xlabel("Predicted risk, EUR-fitted model")
    ax[1].set_ylabel("Observed case fraction")
    ax[1].set_title("B   Calibration", loc="left", fontweight="bold", pad=6)
    leg = ax[1].legend(frameon=False, loc="lower right", handlelength=.9, fontsize=6.4,
                       labelspacing=.25, borderpad=.1, title="O/E", title_fontsize=6.6)
    leg._legend_box.align = "left"
    ax[1].set_xlim(.05, .62); ax[1].set_ylim(.05, .70)

    # (C) flag rate
    xs = np.arange(len(POPS))
    ax[2].bar(xs, [FLAG[p_] for p_ in POPS], color=[PC[p_] for p_ in POPS], width=.62, zorder=3)
    ax[2].axhline(10, color=CB["black"], lw=.9)
    ax[2].annotate("pooled top 10%", (4.45, 10.3), fontsize=6.4, ha="right")
    for i, p_ in enumerate(POPS):
        ax[2].annotate(f"{FLAG[p_]:.1f}", (i, FLAG[p_] + .3), ha="center", fontsize=7, fontweight="bold")
        ax[2].annotate(f"{R2[p_]:.2f}", (i, 0.45), ha="center", fontsize=6.3, color="white",
                       fontweight="bold")
    ax[2].set_xticks(xs); ax[2].set_xticklabels(POPS)
    ax[2].set_ylabel("Flagged at PRS > 1.282\n(% of group)")
    ax[2].set_xlabel("Ancestry (mean tag-causal r$^2$ in bar)")
    ax[2].set_title("C   Who gets flagged", loc="left", fontweight="bold", pad=6)
    ax[2].set_ylim(0, 17)
    save(fig, "figure2_polygenic_risk")

# ---------------------------------------------------------------- Figure 3
def figure3():
    N = pd.read_csv(f"{M7}/nodes.csv"); E = pd.read_csv(f"{M7}/edges.csv")
    PR = {"KRAS": .0377, "BRAF": .0231, "EGFR": .0231, "TYMS": .0231, "CA2": .0158, "MYC": .0158,
          "JAK2": .0158, "TCF7L2": .0158, "BMP4": .0158, "SMAD7": .0158, "TOP1": .0158,
          "PTGS2": .0158, "VEGFA": .0158, "IL12B": .0158, "DES": .0086, "MYH11": .0086,
          "IL23R": .0086, "CR2": .0086}
    kindcol = N.columns[1] if "kind" not in N.columns else "kind"
    kind = dict(zip(N[N.columns[0]], N[kindcol]))
    KC = {"gene": CB["blue"], "disease": CB["red"], "drug": CB["green"],
          "variant": CB["orange"], "pathway": CB["sky"], "phenotype": CB["purple"]}
    ECONF = {"high": (1.0, "-"), "medium": (0.55, "-"), "low": (0.3, ":")}

    G = nx.DiGraph()
    for _, r in N.iterrows(): G.add_node(r[N.columns[0]], kind=r[kindcol])
    for _, r in E.iterrows(): G.add_edge(r.src, r.dst, confidence=r.confidence)

    fig = plt.figure(figsize=(7.4, 4.7))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.52, 1], wspace=0.06, top=0.90)
    axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1])

    pos = nx.spring_layout(G.to_undirected(), seed=7, k=1.18, iterations=450)
    for u, v, d in G.edges(data=True):
        a, ls = ECONF.get(d["confidence"], (.4, "-"))
        axA.annotate("", xy=pos[v], xytext=pos[u],
                     arrowprops=dict(arrowstyle="-", color=CB["grey"], alpha=a, lw=.75, ls=ls))
    for n, d in G.nodes(data=True):
        s = 16 + 5600 * PR.get(n, 0.004)
        axA.scatter(*pos[n], s=s, color=KC.get(d["kind"], CB["grey"]),
                    edgecolor="white", linewidth=.5, zorder=3)
    LABEL = ["KRAS", "BRAF", "EGFR", "TYMS", "TOP1", "MYC", "SMAD7", "TCF7L2", "BMP4",
             "CA2", "CR2", "DES", "MYH11", "IL12B", "Colon adenocarcinoma", "rs6983267"]
    for n in LABEL:
        if n in pos:
            t = "COAD" if n == "Colon adenocarcinoma" else n
            axA.annotate(t, pos[n], xytext=(0, 7.5), textcoords="offset points", ha="center",
                         fontsize=6.6, fontweight="bold" if n in ("KRAS", "MYC", "CA2") else "normal",
                         zorder=5)
    axA.set_axis_off()
    axA.set_title("A   The graph: 54 nodes, 61 edges", loc="left", fontweight="bold", pad=2)
    axA.legend(handles=[Patch(color=c, label=k) for k, c in KC.items()] +
                       [Line2D([], [], color=CB["grey"], ls="-", label="high confidence (25)"),
                        Line2D([], [], color=CB["grey"], ls="-", alpha=.55, label="medium (4)"),
                        Line2D([], [], color=CB["grey"], ls=":", alpha=.5, label="low (32)")],
               loc="lower left", frameon=False, fontsize=6.3, ncol=2, handlelength=1.2,
               bbox_to_anchor=(-0.02, -0.06))
    axA.annotate("node size = PageRank\n(KRAS 0.0377 to CR2 0.0086)", (0.985, 0.985),
                 xycoords="axes fraction", ha="right", va="top", fontsize=6.3, color=CB["grey"])

    # (B) the rs6983267 case
    axB.set_xlim(0, 1); axB.set_ylim(0, 1); axB.set_axis_off()
    axB.set_title("B   One SNP, two gene assignments", loc="left", fontweight="bold", pad=2)
    axB.scatter([.5], [.80], s=300, color=CB["orange"], edgecolor="white", zorder=3)
    axB.annotate("rs6983267", (.5, .80), ha="center", va="center", fontsize=7, fontweight="bold", zorder=4)
    axB.annotate("8q24 locus  |  PRS weight 0.073 (Lab 6)\nGWAS p = 2 x 10$^{-122}$, strongest in the graph",
                 (.5, .99), ha="center", va="top", fontsize=6.6)
    axB.add_patch(FancyArrowPatch((.42, .735), (.22, .50), arrowstyle="-|>", color=CB["red"],
                                  lw=1.3, ls=":", mutation_scale=9))
    axB.add_patch(FancyArrowPatch((.58, .735), (.78, .50), arrowstyle="-|>", color=CB["grey"],
                                  lw=1.0, mutation_scale=9))
    axB.scatter([.22], [.44], s=300, color=CB["blue"], edgecolor="white", zorder=3)
    axB.annotate("MYC", (.22, .44), ha="center", va="center", fontsize=7, color="white",
                 fontweight="bold", zorder=4)
    axB.annotate("what this graph asserts\nOCCURS_IN, confidence LOW\nPageRank 0.0158 (mid-table)",
                 (.22, .355), ha="center", va="top", fontsize=6.5, color=CB["red"])
    axB.annotate("what the GWAS Catalog maps\n$\\it{CASC8}$, $\\it{CCAT2}$, $\\it{PCAT1}$, $\\it{POU5F1B}$\nnone of them $\\it{MYC}$",
                 (.78, .455), ha="center", va="top", fontsize=6.5, color=CB["black"])
    axB.annotate("Same pattern on rs4444235: the graph says $\\it{BMP4}$,\nthe Catalog maps $\\it{DDHD1}$-$\\it{DT}$ and $\\it{MIR5580}$.\n\n"
                 "Lab 6 carries no gene column at all, so a PRS\nscores this variant without naming a gene.\nA graph node has to be named.",
                 (.5, .145), ha="center", va="center", fontsize=6.5,
                 bbox=dict(boxstyle="round,pad=0.45", facecolor=CB["yellow"], alpha=.22, edgecolor="none"))
    save(fig, "figure3_knowledge_graph")

# ---------------------------------------------------------------- Figure 4
def figure4():
    # notebook cell 23 (single cell) and cells 30/39 (spatial)
    ct = [("Epithelial/Tumor", 8678, .633, .364), ("NK", 50, .258, .100), ("DC", 36, .217, .139),
          ("Other", 2976, .192, .171), ("B", 529, .134, .072), ("Fibroblast", 1072, .112, .096),
          ("T", 2354, .056, .025), ("Macrophage", 1189, .048, .030),
          ("Endothelial", 1431, .042, .038), ("Monocyte", 73, .000, .000)]
    dom = [(0, "stroma / immune", 1.22, .80), (1, "stroma / immune", 0.98, .84),
           (2, "tumour nest", 1.52, .99), (3, "tumour nest", 1.31, .96), (4, "tumour nest", 1.65, 1.00),
           (5, "tumour nest", 1.44, .98), (6, "tumour nest", 1.38, .97), (7, "tumour nest", 1.58, .99),
           (8, "tumour nest", 1.41, .98)]

    fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.3), gridspec_kw={"width_ratios": [1.25, 1], "wspace": .42})

    names = [c[0] for c in ct][::-1]
    means = [c[2] for c in ct][::-1]
    fracs = [c[3] for c in ct][::-1]
    ns    = [c[1] for c in ct][::-1]
    ys = np.arange(len(names))
    cols = [CB["red"] if n == "Epithelial/Tumor" else CB["blue"] for n in names]
    ax[0].barh(ys, means, color=cols, height=.6, zorder=3)
    for yy, m, f, nn in zip(ys, means, fracs, ns):
        ax[0].annotate(f"{m:.3f}   {f*100:.1f}% expressing   n={nn:,}", (m + .012, yy),
                       va="center", fontsize=6.4)
    ax[0].axvline(0.26, color=CB["grey"], lw=.8, ls="--")
    ax[0].annotate("0.26", (0.26, 9.55), fontsize=6.3, color=CB["grey"], ha="center", va="bottom")
    ax[0].set_yticks(ys); ax[0].set_yticklabels(names)
    ax[0].set_xlabel("Mean ERBB2 log-normalised expression")
    ax[0].set_xlim(0, 1.03); ax[0].set_ylim(-0.65, 10.2)
    ax[0].set_title("A  $\\it{ERBB2}$ by cell type, 18,388 cells, 10 donors",
                    loc="left", fontweight="bold")

    nest = [d for d in dom if d[1] == "tumour nest"]; stro = [d for d in dom if d[1] != "tumour nest"]
    ax[1].scatter([d[3] * 100 for d in nest], [d[2] for d in nest], s=46, color=CB["red"],
                  edgecolor="white", linewidth=.5, zorder=3, label="tumour nests (domains 2-8)")
    ax[1].scatter([d[3] * 100 for d in stro], [d[2] for d in stro], s=46, color=CB["blue"],
                  edgecolor="white", linewidth=.5, zorder=3, label="stroma / immune (domains 0-1)")
    for d in dom:
        ax[1].annotate(str(d[0]), (d[3] * 100, d[2]), xytext=(0, 7), textcoords="offset points",
                       ha="center", fontsize=6.3)
    ax[1].set_xlabel("Spots expressing ERBB2 (%)")
    ax[1].set_ylabel("Mean ERBB2 expression in domain")
    ax[1].set_xlim(74, 105); ax[1].set_ylim(.88, 1.80)
    ax[1].legend(frameon=False, loc="lower right", fontsize=6.5, handlelength=1.0,
                 bbox_to_anchor=(1.0, -.02))
    ax[1].set_title("B  Spatial domains, Visium, 3,798 spots", loc="left", fontweight="bold")

    save(fig, "figure4_cell_type_resolution")

if __name__ == "__main__":
    print("Figure 1"); figure1()
    print("Figure 2"); figure2()
    print("Figure 3"); figure3()
    print("Figure 4"); figure4()
    print("done ->", OUT)
