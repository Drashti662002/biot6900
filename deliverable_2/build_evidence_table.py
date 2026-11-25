"""Build deliverable_2/evidence_table.csv from existing lab outputs. No new analysis:
every value is read from a file or from a saved notebook output. NA where absent."""
import pandas as pd, numpy as np, os

ROOT = "/Users/drashtibhanushali/biot6900"
T = pd.read_csv(f"{ROOT}/modules/module 7/data/targets_coad.csv").set_index("gene")
D1 = pd.read_csv(f"{ROOT}/deliverable_1/evidence_table.csv").set_index("gene")
SCORE = pd.read_csv(f"{ROOT}/modules/module 6/part3_crc/data/score.csv").set_index("variant")
E = pd.read_csv(f"{ROOT}/modules/module 7/edges.csv")
GL = pd.read_csv(f"{ROOT}/modules/module 7/data/gwas_lookup.csv", dtype={"gwas_p": str})

# PageRank + communities: saved outputs of the committed Lab 7 notebook (bf72daa), cells 14 and 16
PR = {"KRAS":0.0377,"BRAF":0.0231,"EGFR":0.0231,"TYMS":0.0231,"CA2":0.0158,"MYC":0.0158,
      "JAK2":0.0158,"TCF7L2":0.0158,"BMP4":0.0158,"SMAD7":0.0158,"TOP1":0.0158,"PTGS2":0.0158,
      "VEGFA":0.0158,"IL12B":0.0158,"DES":0.0086,"MYH11":0.0086,"IL23R":0.0086,"CR2":0.0086}
COMM = {0:["IL12B","IL23R","JAK2"], 1:["CR2","DES","MYH11"], 2:["BRAF","EGFR"],
        3:["TOP1","TYMS"], 4:["KRAS"], 5:["MYC","TCF7L2"], 6:["BMP4","SMAD7"],
        7:["VEGFA"], 8:["PTGS2"], 9:["CA2"]}
COMM_LABEL = {0:"0 (UC immune: IL12B/IL23R/JAK2 + ustekinumab/tofacitinib)",
  1:"1 (smooth muscle + complement + COAD node: CR2/DES/MYH11)",
  2:"2 (EGFR-RAS-MAPK: BRAF/EGFR + cetuximab/panitumumab/encorafenib)",
  3:"3 (nucleotide synthesis: TYMS/TOP1 + 5-FU/capecitabine/irinotecan)",
  4:"4 (KRAS + G12C/G12D + sotorasib/adagrasib)", 5:"5 (WNT: MYC/TCF7L2 + rs6983267/rs11196172)",
  6:"6 (TGF-beta/BMP: BMP4/SMAD7 + rs4444235/rs4939827)", 7:"7 (VEGFA + bevacizumab)",
  8:"8 (PTGS2 + aspirin)", 9:"9 (CA2 + acetazolamide, carbon dioxide hydration)"}
gene_comm = {g:COMM_LABEL[k] for k,v in COMM.items() for g in v}

PRS_SNP = {"MYC":"rs6983267","SMAD7":"rs4939827","BMP4":"rs4444235","TCF7L2":"rs11196172"}

# Lab 5: per-gene labelled-variant counts recomputed from data/*.tsv (see inventory 1.2)
PGX = {
 "DPYD": ("fluoropyrimidine (5-FU, capecitabine) toxicity; CPIC activity score",
   "105 P/LP in DPYD.tsv, 39 with CADD; *2A rs3918290 CADD 23.7 (no REVEL/AlphaMissense); "
   "*13 rs55886062 CADD 31.0 REVEL 0.92 AM 0.850; c.2846A>T rs67376798 CADD 29.9 REVEL 0.962 AM 0.293 "
   "(excluded, ClinVar Conflicting); HapB3 tag rs56038477 CADD 14.03 counted benign; "
   "causal intronic rs75017182 CADD 0.378"),
 "CYP2C19": ("clopidogrel (not a COAD drug); included as the predictor-scope case",
   "80 ClinVar variants, CADD 68 / REVEL 54 / AlphaMissense 54; *2 rs4244285 CADD 5.686 band 0-10 "
   "(LR+ 0.07, no REVEL/AlphaMissense, synonymous cryptic splice donor); *17 rs12248560 CADD 3.031 band 0-10; "
   "*7 rs72558186 no tool score at all"),
}

def fmt(x, nd=None):
    if x is None or (isinstance(x,float) and np.isnan(x)): return "NA"
    return f"{x:.{nd}f}" if nd is not None else str(x)

rows = []
for g in ["MYC","SMAD7","TCF7L2","BMP4","KRAS","BRAF","TYMS","TOP1","EGFR","VEGFA","PTGS2",
          "CA2","CR2","DES","MYH11","DPYD","CYP2C19"]:
    t = T.loc[g] if g in T.index else None
    in_kg = g in PR
    # axis
    ax = []
    if g in PRS_SNP: ax.append("germline")
    if t is not None and "somatic" in str(t.get("evidence","")): ax.append("somatic")
    if g in PGX or (t is not None and "pharmacology" in str(t.get("evidence",""))): ax.append("pharmacologic")
    if t is not None and "multi-omics" in str(t.get("evidence","")): ax.append("expression")
    axis = "; ".join(ax) if ax else "NA"

    # Part I rank (Deliverable 1 15-gene triage list)
    if g in D1.index: p1 = str(int(D1.loc[g,"lab2_rank"]))
    elif g == "MYH11": p1 = "NA (Lab 2 rank 24; not on the 15-gene triage list)"
    else: p1 = "NA"

    # KG evidence type + confidence, from edges.csv (the gene's ASSOCIATED_WITH edge)
    if in_kg:
        ed = E[(E.src==g) & (E.rel=="ASSOCIATED_WITH")]
        kg_ev = str(t["evidence"]) if t is not None else (ed.source.iloc[0] if len(ed) else "NA")
        kg_conf = ed.confidence.iloc[0] if len(ed) else "NA"
        drugs = E[(E.rel=="TARGETS") & (E.dst==g)].src.tolist()
        if drugs: kg_ev += f"; drug edges: {', '.join(drugs)} (ChEMBL, high)"
    else:
        kg_ev = kg_conf = "NA"

    # PRS
    snp = PRS_SNP.get(g)
    if snp:
        w = SCORE.loc[snp,"weight"]
        in_prs, prs_snp, prs_w = "Y", snp, f"{w:.3f}"
    else:
        in_prs, prs_snp, prs_w = "N", "NA", "NA"

    # read the p-value as written in gwas_lookup.csv, avoiding float round-trip artifacts
    gl = GL[GL.gene==g]
    gp = gl.gwas_p.astype(str).iloc[0] if len(gl) else "NA"
    mapped = "NA"
    if t is not None and isinstance(t.get("gwas_mapped_genes"), str):
        mg = t["gwas_mapped_genes"]
        if g not in mg.split(";"): mapped = mg          # blank where they agree

    pgx_rel, var_scores = PGX.get(g, ("NA","NA"))

    # cell-type evidence
    if g in ("CR2","DES","MYH11"):
        ct = ("inferred from KG community 1 (smooth muscle + complement), no single-cell data; "
              "DES/MYH11 smooth-muscle proteins, CR2 B-cell/FDC receptor")
    elif g == "CA2":
        ct = ("NA (CA2 is in KG community 9 with acetazolamide, NOT in the community-1 "
              "tissue-composition cluster; no single-cell data)")
    else:
        ct = "NA"

    # disease match
    if g in PGX: dm = "COAD drugs (DPYD)" if g=="DPYD" else "no (clopidogrel; CYP2C19 is not a COAD drug gene)"
    else: dm = "COAD (TCGA/CPTAC, GWAS Catalog colorectal cancer)"

    # layers supporting: Part II layers only (germline risk, PGx, network, cell type)
    L = []
    if in_prs=="Y": L.append("germline risk (Lab 6)")
    if g in PGX: L.append("PGx (Lab 5)")
    if in_kg: L.append("network (Lab 7)")
    if g in ("CR2","DES","MYH11"): L.append("cell type (inferred, Lab 7 community)")
    layers = f"{len(L)}: " + "; ".join(L) if L else "0"

    notes = []
    if t is not None:
        notes.append(f"targets_coad.csv: rna_lfc {t['rna_lfc']:.3f}" +
                     (f", prot_lfc {t['prot_lfc']:.3f}" if not pd.isna(t['prot_lfc']) else ", prot_lfc NA") +
                     f", druggability {t['druggability']:.3f}")
    if g in ("CA2","CR2"):
        d1d = D1.loc[g,"druggability_score"]
        notes.append(f"D1 quotes druggability {d1d:.3f} vs Lab 7 refit {t['druggability']:.6f} under the same label")
    if g=="MYC": notes.append("gwas_lookup.csv row 1: PMID 36539618, GCST90129505, Fernandez-Rozadilla C")
    if not in_kg: notes.append("absent from the Lab 7 graph: no node, no edges")
    if g in PRS_SNP: notes.append("score.csv has no gene column; the SNP-to-gene link is Lab 7's")

    rows.append(dict(gene=g, axis=axis, part_I_rank=p1, in_KG="Y" if in_kg else "N",
        KG_PageRank=fmt(PR.get(g),4), KG_community=gene_comm.get(g,"NA"),
        KG_evidence_type=kg_ev, KG_edge_confidence=kg_conf,
        in_PRS=in_prs, PRS_SNP=prs_snp, PRS_weight=prs_w, GWAS_p=gp, GWAS_mapped_genes=mapped,
        PGx_relevance=pgx_rel, variant_scores=var_scores, cell_type_evidence=ct,
        disease_match=dm, layers_supporting=layers, provenance_notes="; ".join(notes) or "NA"))

df = pd.DataFrame(rows)
out = f"{ROOT}/deliverable_2/evidence_table.csv"
df.to_csv(out, index=False)
print("wrote", out, df.shape)
print(df[["gene","axis","part_I_rank","in_KG","KG_PageRank","in_PRS","PRS_weight","GWAS_p","layers_supporting"]].to_string(index=False))
