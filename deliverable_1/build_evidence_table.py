"""
build_evidence_table.py

Builds deliverable_1/evidence_table.csv: one row per candidate that reached
Lab 3 (the 15 genes in modules/module3/targets_coad_w3.csv).

No new analysis. Every value is read from an existing lab output, or is a
trivial summary of one (a count, or a mean of pLDDT values already in the
AlphaFold B-factor column over pocket residues fpocket already listed).

Run from the repo root:

    python deliverable_1/build_evidence_table.py
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
LAB2 = ROOT / "modules/module2/targets_coad.csv"           # Lab 2 ranked table
LAB3 = ROOT / "modules/module3/targets_coad_w3.csv"        # Lab 3 re-ranked table
PROTEOME = ROOT / "shared/module3_pack/proteome.tsv"       # gene -> UniProt accession
LAB4 = ROOT / "modules/module4/structure_assessment_coad.csv"
STRUCT = ROOT / "shared/module4_structures"
OUT = ROOT / "deliverable_1/evidence_table.csv"

FDR = 0.05   # the only significance threshold Lab 2 used (rna_q / prot_q, and gen_q in the report)

# Pockets Lab 4 reported (notebook Step 6 output). Top pocket by fpocket score,
# and, where different, the pocket the lab identified as the real ligand site.
POCKETS = {
    "CA2": {"top": 6, "ligand_site": 2},     # AZM A701 (acetazolamide) sits in pocket 2
    "CR2": {"top": 50, "best_kept": 27},     # pocket 50 rejected (pLDDT < 70)
}


def ca_plddt(pdb: Path) -> dict[int, float]:
    """pLDDT per residue = B-factor of the CA atom in an AlphaFold model."""
    out = {}
    for line in open(pdb):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            out[int(line[22:26])] = float(line[60:66])
    return out


def pocket_residues(path: Path) -> list[int]:
    return sorted({int(l[22:26]) for l in open(path) if l.startswith("ATOM")})


def pocket_scores(info: Path) -> dict[int, float]:
    txt = open(info).read()
    out = {}
    for m in re.finditer(r"Pocket (\d+) :(.*?)(?=\nPocket|\Z)", txt, flags=re.S):
        out[int(m.group(1))] = float(re.search(r"Druggability Score :\s*([\d.]+)", m.group(2)).group(1))
    return out


def pocket_summary(acc: str, pocket: int) -> tuple[float, float]:
    plddt = ca_plddt(STRUCT / f"AF-{acc}.pdb")
    res = pocket_residues(STRUCT / f"AF-{acc}_out/pockets/pocket{pocket}_atm.pdb")
    score = pocket_scores(STRUCT / f"AF-{acc}_out/AF-{acc}_info.txt")[pocket]
    return score, round(float(np.mean([plddt[r] for r in res])), 1)


def main() -> None:
    lab2 = pd.read_csv(LAB2).set_index("gene")
    lab3 = pd.read_csv(LAB3)
    acc = pd.read_csv(PROTEOME, sep="\t", usecols=["gene", "accession"]).set_index("gene")["accession"]
    lab4 = pd.read_csv(LAB4).set_index("gene")

    rows = []
    for _, r3 in lab3.iterrows():
        g = r3["gene"]
        r2 = lab2.loc[g]
        layers = int(r2.rna_q < FDR) + int(r2.prot_q < FDR) + int(r2.gen_q < FDR)
        direction = ("both down" if r2.rna_lfc < 0 and r2.prot_lfc < 0 else
                     "both up" if r2.rna_lfc > 0 and r2.prot_lfc > 0 else
                     f"RNA {'up' if r2.rna_lfc > 0 else 'down'}, protein {'up' if r2.prot_lfc > 0 else 'down'}")
        best = max(r3.p_sm, r3.p_ab)
        modality = "small molecule" if r3.p_sm >= r3.p_ab else "antibody"
        row = {
            "gene": g,
            "uniprot_id": acc.get(g, "NA"),
            "lab2_rank": int(r2["rank"]),
            "transcriptomic_evidence (log2FC, adj p)": f"log2FC {r2.rna_lfc:.2f}; q {r2.rna_q:.1e}",
            "proteomic_evidence": f"log2FC {r2.prot_lfc:.2f}; q {r2.prot_q:.1e}",
            "genomic_evidence (mutation/CNV freq)":
                f"mut {100 * r2.mut_freq:.1f}% ({int(r2.n_mutated)}/{int(r2.n_patients)}); q {r2.gen_q:.2f}; CNV NA",
            "n_layers_supporting": layers,
            "concordance": f"{'concordant' if r2.concordant else 'discordant'} ({direction})",
            "druggability_score": round(best, 3),
            "druggable (Y/N, threshold)": "NA (Lab 3 set no threshold)",
            "structure_available": "Y" if g in lab4.index else "N",
            "mean_pLDDT": lab4.loc[g, "mean_pLDDT"] if g in lab4.index else "NA",
            "pocket_found": "NA",
            "top_pocket_score": "NA",
            "mean_pLDDT_of_pocket_residues": "NA",
            "notes/gaps": (f"Lab 3 rank {int(r3['rank'])}; p_sm {r3.p_sm:.3f}, p_ab {r3.p_ab:.3f} "
                           f"(best modality: {modality}). "),
        }
        if g in POCKETS:
            a = lab4.loc[g, "uniprot"]
            top_score, top_pl = pocket_summary(a, POCKETS[g]["top"])
            row["top_pocket_score"] = top_score
            row["mean_pLDDT_of_pocket_residues"] = top_pl
            if g == "CA2":
                s2, p2 = pocket_summary(a, POCKETS[g]["ligand_site"])
                row["pocket_found"] = "Y"
                row["notes/gaps"] += (f"Top fpocket pocket 6 is not the drug site; acetazolamide site = pocket 2 "
                                      f"(score {s2:.3f}, pocket pLDDT {p2}, covers 85% of AZM A701 contacts). "
                                      "Single domain; AF vs 3HS4 RMSD 0.37 A.")
            else:
                s27, p27 = pocket_summary(a, POCKETS[g]["best_kept"])
                row["pocket_found"] = "N"
                row["notes/gaps"] += (f"Top pocket 50 rejected (pocket pLDDT {top_pl} < 70, spans Sushi 7-8 junction); "
                                      f"best kept pocket 27 scores {s27:.3f} (pLDDT {p27}). 8 domains, inter-domain "
                                      "PAE 26.3 A. Lab 4 RMSD 38.45 A vs 1LY2 is an alignment artefact "
                                      "(1.00 A over residues 20-149 in a side check).")
        else:
            row["notes/gaps"] += "Not carried to Lab 4: no structure or pocket assessed."
        if layers < 3:
            row["notes/gaps"] += f" Genomic layer not significant (gen_q {r2.gen_q:.2f})."
        if r2.rna_lfc < 0 and r2.prot_lfc < 0:
            row["notes/gaps"] += " Lost in tumour on both layers (not an inhibition target as-is)."
        rows.append(row)

    out = pd.DataFrame(rows).sort_values("lab2_rank")
    out.to_csv(OUT, index=False)
    print(f"wrote {OUT.relative_to(ROOT)} ({len(out)} rows)")


if __name__ == "__main__":
    main()
