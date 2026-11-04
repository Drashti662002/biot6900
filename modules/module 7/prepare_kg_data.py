"""Build data/targets_coad.csv for the Module 7 knowledge graph.

One row per gene node, with the columns the starter notebook reads
(gene, ensembl, gwas_p, rna_lfc, prot_lfc, druggability) plus provenance.

Sources, all from earlier modules except the GWAS lookup:
  - rna_lfc / prot_lfc / mutation p: Module 2 per-layer COAD tables
    (TCGA-COAD RNA, CPTAC-2 protein, CPTAC-2 somatic mutation recurrence).
  - ensembl, druggability: Module 3 data pack. druggability = max(p_sm, p_ab)
    from the Module 3 tractability prior (ESM C embeddings + logistic regression,
    same model and settings), refit with every gene scored here held out.
  - gwas_p: GWAS Catalog REST API v2, best p-value for the lead SNP's
    association with the named trait. Cached in data/gwas_lookup.csv.

Run from this folder:  python prepare_kg_data.py
"""
import json
import re
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).resolve().parent
SHARED = HERE.parent.parent / "shared"
OMICS, PACK = SHARED / "module2_omics", SHARED / "module3_pack"
OUT = HERE / "data"

# gene -> why it is in the graph (becomes the ASSOCIATED_WITH edge's `source`)
GENES = {
    # Module 2 top-ranked, RNA+protein concordant, all lower in tumour
    "CA2": "multi-omics (TCGA/CPTAC)", "DES": "multi-omics (TCGA/CPTAC)",
    "CR2": "multi-omics (TCGA/CPTAC)", "MYH11": "multi-omics (TCGA/CPTAC)",
    # recurrently mutated drivers with approved CRC drugs
    "KRAS": "somatic mutation (CPTAC)", "BRAF": "somatic mutation (CPTAC)",
    # approved CRC drug targets
    "EGFR": "pharmacology", "VEGFA": "pharmacology", "TYMS": "pharmacology",
    "TOP1": "pharmacology", "PTGS2": "pharmacology",
    # germline CRC GWAS loci from the Module 6 PRS (PGS000148)
    "MYC": "GWAS", "SMAD7": "GWAS", "BMP4": "GWAS", "TCF7L2": "GWAS",
}
RELATED_GENES = ["IL23R", "JAK2", "IL12B"]   # ulcerative colitis GWAS genes

# (gene, lead SNP, trait regex). Gene is the literature candidate gene at the
# locus; the GWAS Catalog's mapped gene is recorded alongside for transparency.
GWAS_SNPS = [
    ("MYC", "rs6983267", "colorectal cancer"),
    ("SMAD7", "rs4939827", "colorectal cancer"),
    ("BMP4", "rs4444235", "colorectal cancer"),
    ("TCF7L2", "rs11196172", "colorectal cancer"),
    ("IL23R", "rs11209026", "ulcerative colitis"),
    ("JAK2", "rs10758669", "ulcerative colitis"),
    ("IL12B", "rs6871626", "ulcerative colitis"),
]


def gwas_best(rs, trait):
    rows, page = [], 0
    while True:
        url = (f"https://www.ebi.ac.uk/gwas/rest/api/v2/associations"
               f"?rs_id={rs}&size=200&page={page}")
        d = json.load(urllib.request.urlopen(url, timeout=60))
        rows += d.get("_embedded", {}).get("associations", [])
        if page + 1 >= d["page"]["totalPages"]:
            break
        page += 1
    # every trait must match, so pooled cross-disease studies are not counted
    hits = [a for a in rows
            if all(re.fullmatch(trait, t["efo_trait"], re.I) for t in a["efo_traits"])]
    best = min(hits, key=lambda a: a["p_value"])
    return dict(gwas_p=best["p_value"], gwas_pmid=best["pubmed_id"],
                gwas_study=best["accession_id"], gwas_author=best["first_author"],
                gwas_mapped_genes=";".join(best["mapped_genes"]))


def gwas_table():
    cache = OUT / "gwas_lookup.csv"
    if cache.exists():
        return pd.read_csv(cache)
    rows = [dict(gene=g, gwas_snp=rs, gwas_trait=t, **gwas_best(rs, t)) for g, rs, t in GWAS_SNPS]
    df = pd.DataFrame(rows)
    df.to_csv(cache, index=False)
    return df


def druggability(genes):
    # Module 3 prior: StandardScaler + LogisticRegression(C=0.05) on ESM C 300M
    # embeddings, one model per modality. Every gene scored here is removed
    # from training so none is scored by a model that saw its own label.
    panel = pd.read_csv(PACK / "panel.tsv", sep="\t")
    X_panel = np.load(PACK / "panel_emb_main.npy").astype("float32")
    proteome = pd.read_csv(PACK / "proteome.tsv", sep="\t", usecols=["gene", "accession", "ensembl_gene"])
    X_prot = np.load(PACK / "proteome_emb_main.npy").astype("float32")

    keep = ~panel["gene"].isin(genes).to_numpy()
    model = lambda: make_pipeline(StandardScaler(), LogisticRegression(C=0.05, max_iter=5000))
    sm = model().fit(X_panel[keep], panel.loc[keep, "sm_clinical"].astype(int))
    ab = model().fit(X_panel[keep], panel.loc[keep, "ab_clinical"].astype(int))

    row = {g: i for i, g in enumerate(proteome["gene"])}
    X = X_prot[[row[g] for g in genes]]
    return pd.DataFrame({"gene": genes,
                         "ensembl": proteome.set_index("gene").loc[genes, "ensembl_gene"].values,
                         "uniprot": proteome.set_index("gene").loc[genes, "accession"].values,
                         "p_sm": sm.predict_proba(X)[:, 1], "p_ab": ab.predict_proba(X)[:, 1]})


def main():
    OUT.mkdir(exist_ok=True)
    all_genes = list(GENES) + RELATED_GENES

    rna = pd.read_csv(OMICS / "coad_transcriptomics.tsv", sep="\t").set_index("gene")
    prot = pd.read_csv(OMICS / "coad_proteomics.tsv", sep="\t").set_index("gene")
    mut = pd.read_csv(OMICS / "coad_genomics.tsv", sep="\t").set_index("gene")

    df = druggability(all_genes)
    df["druggability"] = df[["p_sm", "p_ab"]].max(axis=1)
    df = df.merge(gwas_table(), on="gene", how="left")
    df["evidence"] = df["gene"].map(GENES).fillna("GWAS (ulcerative colitis)")
    df["disease"] = np.where(df["gene"].isin(RELATED_GENES), "Ulcerative colitis", "Colon adenocarcinoma")

    for col, tab, src in [("rna", rna, "log2fc"), ("prot", prot, "log2fc")]:
        df[f"{col}_lfc"] = df["gene"].map(tab[src])
        df[f"{col}_q"] = df["gene"].map(tab["qval"])
    df["mut_n"] = df["gene"].map(mut["n_mutated"])
    df["mut_q"] = df["gene"].map(mut["qval"])
    # COAD expression/mutation evidence does not apply to the related disease
    df.loc[df["disease"] != "Colon adenocarcinoma", ["rna_lfc", "rna_q", "prot_lfc", "prot_q", "mut_n", "mut_q"]] = np.nan

    cols = ["gene", "ensembl", "gwas_p", "rna_lfc", "prot_lfc", "druggability",
            "disease", "evidence", "uniprot", "p_sm", "p_ab", "rna_q", "prot_q", "mut_n", "mut_q",
            "gwas_snp", "gwas_trait", "gwas_pmid", "gwas_study", "gwas_author", "gwas_mapped_genes"]
    print("\nrelated-disease genes (entered by hand in the notebook's RELATED):\n",
          df.loc[df["disease"] != "Colon adenocarcinoma", ["gene", "ensembl", "gwas_p", "druggability"]].round(3))
    df = df.loc[df["disease"] == "Colon adenocarcinoma", cols]
    df.to_csv(OUT / "targets_coad.csv", index=False)
    print(f"wrote {OUT / 'targets_coad.csv'} ({len(df)} genes)")

    # sanity check: the refit prior against the published Module 3 values
    w3 = pd.read_csv(HERE.parent / "module3" / "targets_coad_w3.csv").set_index("gene")
    chk = df.set_index("gene").join(w3[["p_sm", "p_ab"]], rsuffix="_module3", how="inner")
    print("\nrefit prior vs Module 3:\n", chk[["p_sm", "p_sm_module3", "p_ab", "p_ab_module3"]].round(3))

    # somatic hotspot counts used for the VARIANTS layer
    site = pd.read_csv(OMICS / "raw" / "cptac_coad_mutation_site.cbt", sep="\t", index_col=0)
    for v in ["KRAS_p.G12D", "KRAS_p.G12C", "BRAF_p.V600E"]:
        print(f"{v}: {int(site.loc[v].sum())} of {site.shape[1]} tumours")


if __name__ == "__main__":
    main()
