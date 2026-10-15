"""
prepare_coad_data.py

This builds the three gene-level tables that Part 3 of the notebook reads.

Why I needed it: Part 3 normally gives you three ready-made files from Canvas
(ad_transcriptomics.tsv, ad_proteomics.tsv, ad_gwas.tsv). I switched the disease
to colon adenocarcinoma, which my professor approved, so there was no Canvas
file to download. I had to make the equivalent three tables myself from real
public data. I kept this out of the notebook so the notebook stays quick to run
and doesn't re-download 175 MB every time.

Run it once from the repo root:

    python prepare_coad_data.py

It writes:

    data/coad_transcriptomics.tsv   gene, log2fc, pval (+ a few extra columns)
    data/coad_proteomics.tsv        gene, log2fc, pval (+ a few extra columns)
    data/coad_genomics.tsv          gene, neglog10p (+ a few extra columns)
    data/PROVENANCE.md              where every number came from

I used the same column names the professor's files use (gene / log2fc / pval /
neglog10p) so the Part 3 steps work without changing his code.

It also writes three smaller files that I only use for context in the report:
the matched-sample RNA-protein correlation, the KRAS variant spectrum, and an
independent proteomics cohort for validation.

Where the data comes from (all open access, no data use agreement):

1. RNA - TCGA Colon Adenocarcinoma (TCGA-COAD), TCGA Network, Nature
   2012;487:330-337. STAR log2(TPM+1) from the GDC, mirrored by UCSC Xena.
   471 primary tumours and 41 adjacent normals.

2. Protein - CPTAC-2 Prospective Colon Cancer, Vasaikar et al., Cell
   2019;177:1035-1049.e19 (PMID 31031003). PNNL TMT gene-level log ratios,
   mirrored by LinkedOmics. 96 tumour/normal pairs from the same patients.

3. Mutations - the whole-exome data from that same CPTAC study, so the same
   patients as the protein layer. WUSM GATK gene-level calls, also from
   LinkedOmics. 106 patients.

Gene coding lengths come from the UCSC hg38 refFlat file. I use them to correct
the mutation layer for gene length.
"""

from __future__ import annotations

import glob
import gzip
import os
import sys
import urllib.request
import warnings

import numpy as np
import pandas as pd
from scipy import stats

# ------------------------------------------------------------------- settings

RAW = "data/raw"          # the big downloads go here, and git ignores this folder
OUT = "data"              # the small tables the notebook reads

LINKEDOMICS = "https://www.linkedomics.org/data_download/CPTAC-COAD"
XENA = "https://gdc-hub.s3.us-east-1.amazonaws.com/download"

DOWNLOADS = {
    # what I call it locally -> where it comes from
    "TCGA-COAD.star_tpm.tsv.gz": f"{XENA}/TCGA-COAD.star_tpm.tsv.gz",
    "gencode.v36.probemap": f"{XENA}/gencode.v36.annotation.gtf.gene.probemap",
    "refFlat.txt.gz": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz",
    "cptac_coad_prot_tumor.cct": (
        f"{LINKEDOMICS}/Human__CPTAC_COAD__PNNL__Proteome__TMT__03_01_2017"
        "__BCM__Gene__PNNL_Tumor_TMT_UnsharedLogRatio.cct"
    ),
    "cptac_coad_prot_normal.cct": (
        f"{LINKEDOMICS}/Human__CPTAC_COAD__PNNL__Proteome__TMT__03_01_2017"
        "__BCM__Gene__PNNL_Normal_TMT_UnsharedLogRatio.cct"
    ),
    "cptac_coad_prot_tn_log2fc.cct": (
        f"{LINKEDOMICS}/Human__CPTAC_COAD__PNNL__Proteome__TMT__03_01_2017"
        "__BCM__Gene__Tumor_Normal_log2FC.cct"
    ),
    "cptac_coad_rnaseq_tumor.cct": (
        f"{LINKEDOMICS}/Human__CPTAC_COAD__UNC__RNAseq__HiSeq_RNA__03_01_2017"
        "__BCM__Gene__BCM_RSEM_UpperQuartile_log2.cct.gz"
    ),
    "cptac_coad_mutation_gene.cbt": (
        f"{LINKEDOMICS}/Human__CPTAC_COAD__WUSM__Mutation__GAIIx__03_01_2017"
        "__BCM__Gene__GATK_Pipeline.cbt"
    ),
    "cptac_coad_mutation_site.cbt": (
        f"{LINKEDOMICS}/Human__CPTAC_COAD__WUSM__Mutation__GAIIx__03_01_2017"
        "__BCM__Site__GATK_Pipeline.cbt"
    ),
}

# Below this mean log2(TPM+1) in both groups, a gene is basically not expressed
# and its "fold change" is just noise on tiny numbers, so I drop it.
MIN_EXPRESSION = 1.0

# I didn't want to report a paired p-value off a handful of pairs.
MIN_PAIRS = 20

# I fill this in as each step runs and then paste it into PROVENANCE.md, so the
# documentation can't drift away from what the code actually did.
STATS: dict[str, object] = {"min_pairs": MIN_PAIRS}


def log(msg: str) -> None:
    print(msg, flush=True)


def fetch(name: str, url: str) -> str:
    """Download url into data/raw/name, unless it's already there."""
    path = os.path.join(RAW, name)
    if os.path.exists(path) and os.path.getsize(path) > 0:
        log(f"  [cached] {name} ({os.path.getsize(path)/1e6:.1f} MB)")
        return path
    log(f"  [get]    {name} <- {url}")
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                          "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/140.0 Safari/537.36",
            "Referer": LINKEDOMICS + "/",
        },
    )
    tmp = path + ".part"
    with urllib.request.urlopen(req, timeout=600) as r, open(tmp, "wb") as f:
        while chunk := r.read(1 << 20):
            f.write(chunk)
    # A couple of the LinkedOmics files are served gzipped even though I save
    # them under a plain name, so unzip here and the readers don't have to care.
    if url.endswith(".gz") and not name.endswith(".gz"):
        with gzip.open(tmp, "rb") as gz, open(tmp + ".out", "wb") as out:
            while chunk := gz.read(1 << 20):
                out.write(chunk)
        os.replace(tmp + ".out", tmp)
    os.replace(tmp, path)
    log(f"           done ({os.path.getsize(path)/1e6:.1f} MB)")
    return path


# ------------------------------------------------------------------ RNA layer

def build_transcriptomics() -> pd.DataFrame:
    """TCGA-COAD tumour vs adjacent normal, summarised per gene symbol."""
    log("\n[1/3] TRANSCRIPTOMICS  (TCGA-COAD, tumour vs adjacent normal)")
    tpm_path = fetch("TCGA-COAD.star_tpm.tsv.gz", DOWNLOADS["TCGA-COAD.star_tpm.tsv.gz"])
    map_path = fetch("gencode.v36.probemap", DOWNLOADS["gencode.v36.probemap"])

    # Xena gives these as log2(TPM + 1), with versioned Ensembl IDs as the rows.
    with gzip.open(tpm_path, "rt") as fh:
        expr = pd.read_csv(fh, sep="\t", index_col=0)
    log(f"  matrix: {expr.shape[0]} Ensembl IDs x {expr.shape[1]} samples")

    # The 4th field of a TCGA barcode says what the sample is: 01 is primary
    # tumour, 11 is solid tissue normal. I leave out 02 and 06 (recurrent and
    # metastatic) because they aren't part of a tumour vs normal comparison.
    def sample_type(barcode: str) -> str:
        parts = barcode.split("-")
        return parts[3][:2] if len(parts) > 3 else "??"

    types = pd.Series({c: sample_type(c) for c in expr.columns})
    tumor_cols = list(types.index[types == "01"])
    normal_cols = list(types.index[types == "11"])
    log(f"  sample types: {dict(types.value_counts())}")
    log(f"  using {len(tumor_cols)} primary tumours vs {len(normal_cols)} adjacent normals")
    if len(normal_cols) < 5:
        sys.exit("ERROR: too few adjacent normals to compare against.")

    # Map Ensembl IDs to gene symbols. Some symbols have more than one Ensembl
    # ID, so I keep the highest-expressed one instead of averaging them. If I
    # averaged, a dead probe would drag down a real signal.
    probemap = pd.read_csv(map_path, sep="\t", usecols=["id", "gene"]).dropna()
    expr = expr.loc[expr.index.intersection(probemap["id"])]
    symbol = probemap.set_index("id")["gene"].reindex(expr.index)
    order = expr.mean(axis=1).sort_values(ascending=False).index
    expr = expr.loc[order]
    symbol = symbol.loc[order]
    keep = ~symbol.duplicated(keep="first")
    expr, symbol = expr.loc[keep.values], symbol.loc[keep.values]
    expr.index = symbol.str.strip().str.upper().values
    expr.index.name = "gene"
    log(f"  collapsed to {expr.shape[0]} unique gene symbols")

    tumor = expr[tumor_cols]
    normal = expr[normal_cols]
    tmean, nmean = tumor.mean(axis=1), normal.mean(axis=1)

    expressed = (tmean >= MIN_EXPRESSION) | (nmean >= MIN_EXPRESSION)
    log(f"  {int(expressed.sum())} genes pass the expression cutoff "
        f"(mean log2(TPM+1) >= {MIN_EXPRESSION} in tumour or in normal)")
    tumor, normal = tumor[expressed], normal[expressed]
    tmean, nmean = tmean[expressed], nmean[expressed]

    # The values are already logged, so subtracting the means gives log2FC.
    log2fc = tmean - nmean
    # Mann-Whitney because these are different patients (so not paired) and I
    # don't want to assume the data is normally distributed.
    pval = stats.mannwhitneyu(
        tumor.to_numpy(), normal.to_numpy(), axis=1, alternative="two-sided"
    ).pvalue

    out = pd.DataFrame(
        {
            "gene": tumor.index,
            "log2fc": log2fc.to_numpy(),
            "pval": pval,
            "mean_tumor_log2tpm": tmean.to_numpy(),
            "mean_normal_log2tpm": nmean.to_numpy(),
            "n_tumor": len(tumor_cols),
            "n_normal": len(normal_cols),
        }
    )
    out["qval"] = benjamini_hochberg(out["pval"].to_numpy())
    log(f"  -> {len(out)} genes;  {(out['qval'] < 0.05).sum()} with FDR < 0.05")
    STATS.update(n_tumor=len(tumor_cols), n_normal=len(normal_cols), n_tx=len(out))
    return out.sort_values("gene").reset_index(drop=True)


# -------------------------------------------------------------- protein layer

def build_proteomics() -> pd.DataFrame:
    """CPTAC-2 colon, tumour vs adjacent normal, paired within each patient."""
    log("\n[2/3] PROTEOMICS  (CPTAC-2 colon, matched tumour/normal pairs)")
    t_path = fetch("cptac_coad_prot_tumor.cct", DOWNLOADS["cptac_coad_prot_tumor.cct"])
    n_path = fetch("cptac_coad_prot_normal.cct", DOWNLOADS["cptac_coad_prot_normal.cct"])
    ref_path = fetch("cptac_coad_prot_tn_log2fc.cct", DOWNLOADS["cptac_coad_prot_tn_log2fc.cct"])

    tumor = pd.read_csv(t_path, sep="\t", index_col=0)
    normal = pd.read_csv(n_path, sep="\t", index_col=0)
    log(f"  tumour: {tumor.shape}   normal: {normal.shape}")

    # Normal columns are named like 01CO005N, tumour columns like 01CO005, so
    # dropping the trailing N lines the two files up patient by patient.
    normal.columns = [c[:-1] if c.endswith("N") else c for c in normal.columns]
    pairs = sorted(set(tumor.columns) & set(normal.columns))
    log(f"  {len(pairs)} patients have both a tumour and an adjacent-normal channel")

    genes = tumor.index.intersection(normal.index)
    diff = tumor.loc[genes, pairs].to_numpy() - normal.loc[genes, pairs].to_numpy()
    # diff is now tumour minus normal for each patient, which is the same thing
    # the study's own Tumor_Normal_log2FC file contains. I check that below.

    n_pairs = np.sum(~np.isnan(diff), axis=1)
    with warnings.catch_warnings():
        # Some genes aren't quantified in any pair, which gives an all-NaN row
        # and a warning. Those rows get dropped anyway.
        warnings.simplefilter("ignore", RuntimeWarning)
        log2fc = np.nanmean(diff, axis=1)

    # Wilcoxon signed-rank because the pairing is real (same patient), which
    # makes it the right test and also gives more power than an unpaired one.
    pval = np.full(len(genes), np.nan)
    for i in range(len(genes)):
        row = diff[i][~np.isnan(diff[i])]
        if len(row) >= MIN_PAIRS and np.any(row != 0):
            pval[i] = stats.wilcoxon(row, zero_method="wilcox").pvalue

    out = pd.DataFrame(
        {
            "gene": [str(g).strip().upper() for g in genes],
            "log2fc": log2fc,
            "pval": pval,
            "n_pairs": n_pairs,
        }
    ).dropna(subset=["pval"])
    out["qval"] = benjamini_hochberg(out["pval"].to_numpy())
    log(f"  -> {len(out)} genes with at least {MIN_PAIRS} usable pairs; "
        f"{(out['qval'] < 0.05).sum()} with FDR < 0.05")

    # Check: do I get the same log2FC the study published? If not, I've made a
    # mistake somewhere and shouldn't trust this layer.
    ref = pd.read_csv(ref_path, sep="\t", index_col=0)
    ref_lfc = ref.mean(axis=1)
    shared = out.set_index("gene").index.intersection(
        [str(g).strip().upper() for g in ref_lfc.index]
    )
    ref_lfc.index = [str(g).strip().upper() for g in ref_lfc.index]
    a = out.set_index("gene").loc[shared, "log2fc"]
    b = ref_lfc.loc[shared]
    r = np.corrcoef(a.to_numpy(), b.to_numpy())[0, 1]
    log(f"  CHECK against the study's own Tumor_Normal_log2FC ({len(shared)} shared "
        f"genes): Pearson r = {r:.4f}")
    if r < 0.95:
        log("  !! That's lower than I expected. Something is wrong with this layer.")
    STATS.update(n_pairs=len(pairs), n_pr=len(out), prot_r=r, n_shared=len(shared))
    return out.sort_values("gene").reset_index(drop=True)


# -------------------------------------------------------------- mutation layer

def cds_lengths(refflat_path: str) -> pd.Series:
    """Longest coding sequence length per gene symbol, from UCSC refFlat."""
    cols = ["geneName", "name", "chrom", "strand", "txStart", "txEnd",
            "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds"]
    with gzip.open(refflat_path, "rt") as fh:
        rf = pd.read_csv(fh, sep="\t", names=cols, low_memory=False)
    # Only the main chromosomes. The alt and random contigs repeat genes.
    rf = rf[~rf["chrom"].str.contains("_")]

    best: dict[str, int] = {}
    for gene, cds_s, cds_e, ex_s, ex_e in zip(
        rf["geneName"], rf["cdsStart"], rf["cdsEnd"], rf["exonStarts"], rf["exonEnds"]
    ):
        if cds_e <= cds_s:      # this transcript doesn't code for anything
            continue
        starts = [int(x) for x in str(ex_s).rstrip(",").split(",") if x]
        ends = [int(x) for x in str(ex_e).rstrip(",").split(",") if x]
        # CDS length is the exon bases that fall inside the coding region.
        length = sum(
            max(0, min(e, cds_e) - max(s, cds_s)) for s, e in zip(starts, ends)
        )
        key = str(gene).strip().upper()
        if length > best.get(key, 0):
            best[key] = length
    return pd.Series(best, name="cds_len")


def build_genomics() -> pd.DataFrame:
    """CPTAC-2 colon exome data, turned into a mutation recurrence score."""
    log("\n[3/3] GENOMICS  (CPTAC-2 colon exomes, somatic mutation recurrence)")
    mut_path = fetch("cptac_coad_mutation_gene.cbt", DOWNLOADS["cptac_coad_mutation_gene.cbt"])
    rf_path = fetch("refFlat.txt.gz", DOWNLOADS["refFlat.txt.gz"])

    mut = pd.read_csv(mut_path, sep="\t", index_col=0)
    mut.index = [str(g).strip().upper() for g in mut.index]
    n_patients = mut.shape[1]
    log(f"  gene x patient table of 0s and 1s: {mut.shape} "
        f"(values {sorted(pd.unique(mut.to_numpy().ravel()))})")

    n_mut = (mut > 0).sum(axis=1)
    freq = n_mut / n_patients

    lengths = cds_lengths(rf_path)
    log(f"  got CDS lengths for {len(lengths)} gene symbols (UCSC hg38 refFlat)")

    df = pd.DataFrame({"gene": mut.index, "n_mutated": n_mut.to_numpy(),
                       "mut_freq": freq.to_numpy()})
    df["cds_len"] = df["gene"].map(lengths)
    have = df["cds_len"].notna()
    log(f"  {int(have.sum())}/{len(df)} mutated genes have a CDS length")

    # I started with plain mutation frequency and the top of the list was TTN,
    # MUC16 and OBSCN, which are just very long genes that pick up passenger
    # mutations. So instead I work out one background mutation rate per coding
    # base per patient for the whole cohort, scale it by each gene's CDS length
    # to get how often you'd expect that gene to be hit by chance, and test the
    # observed count against that. It's the same idea as MutSigCV but much
    # simpler. It's a correction for gene length, not a proper driver-gene
    # method, and I don't claim it is one.
    total_events = df.loc[have, "n_mutated"].sum()
    total_bases = df.loc[have, "cds_len"].sum()
    mu = total_events / (total_bases * n_patients)
    log(f"  background rate mu = {mu:.3e} mutated-gene events per coding base per patient")

    # chance that this gene is hit in any one patient, under the background rate
    p_gene = 1.0 - np.exp(-mu * df["cds_len"])
    # one-sided, because I only care about genes mutated more often than expected
    pval = stats.binom.sf(df["n_mutated"] - 1, n_patients, p_gene)
    df["pval"] = pval
    df["neglog10p"] = -np.log10(np.clip(pval, 1e-300, 1.0))
    df.loc[~have, ["pval", "neglog10p"]] = np.nan

    df = df.dropna(subset=["neglog10p"])
    df["qval"] = benjamini_hochberg(df["pval"].to_numpy())
    df["n_patients"] = n_patients
    log(f"  -> {len(df)} genes scored; {(df['qval'] < 0.05).sum()} significant at FDR < 0.05")
    top = df.nlargest(10, "neglog10p")[["gene", "n_mutated", "mut_freq", "cds_len", "neglog10p"]]
    log("  top 10 after the length correction:")
    log(top.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    STATS.update(n_patients=n_patients, mu=mu, n_gw=len(df))
    return df.sort_values("gene").reset_index(drop=True)


# -------------------------------------------------- extra context I use in the report

def build_matched_rna_protein_correlation() -> pd.DataFrame:
    """RNA vs protein correlation across CPTAC tumours, like Part 1 does.

    CPTAC measured RNA and protein on the same colon tumours, so for those
    samples I can ask the Part 1 question: does protein follow RNA across
    patients? My Part 3 join can't ask that, because TCGA and CPTAC are
    different people. Putting a real number on the gap seemed better than just
    saying there is one.
    """
    log("\n[4/6] CONTEXT: RNA vs protein in the same CPTAC tumours")
    rna_path = fetch("cptac_coad_rnaseq_tumor.cct", DOWNLOADS["cptac_coad_rnaseq_tumor.cct"])
    prot_path = os.path.join(RAW, "cptac_coad_prot_tumor.cct")

    rna = pd.read_csv(rna_path, sep="\t", index_col=0)
    prot = pd.read_csv(prot_path, sep="\t", index_col=0)
    rna.index = [str(g).strip().upper() for g in rna.index]
    prot.index = [str(g).strip().upper() for g in prot.index]

    samples = sorted(set(rna.columns) & set(prot.columns))
    genes = sorted(set(rna.index) & set(prot.index))
    log(f"  {len(samples)} tumours measured on both layers; {len(genes)} shared genes")

    R = rna.loc[genes, samples]
    P = prot.loc[genes, samples]
    rows = []
    for g in genes:
        a, b = R.loc[g].to_numpy(float), P.loc[g].to_numpy(float)
        ok = ~(np.isnan(a) | np.isnan(b))
        if ok.sum() < 20 or np.std(a[ok]) == 0 or np.std(b[ok]) == 0:
            continue
        rows.append((g, stats.spearmanr(a[ok], b[ok]).statistic, int(ok.sum())))
    out = pd.DataFrame(rows, columns=["gene", "spearman_r", "n_samples"])
    med = out["spearman_r"].median()
    log(f"  -> {len(out)} genes; median Spearman r = {med:.3f} "
        f"(from {out['spearman_r'].min():.2f} to {out['spearman_r'].max():.2f})")
    STATS.update(n_matched_samples=len(samples), n_corr_genes=len(out), median_corr=med)
    return out


def build_kras_variant_spectrum() -> pd.DataFrame:
    """Which KRAS mutations actually show up in this cohort, and how often.

    KRAS is the gene I set out to look at, and it would have been easy to
    assume this was a G12C cohort just because G12C is the one with drugs
    against it. This table is me checking instead of assuming.
    """
    log("\n[5/6] CONTEXT: which KRAS variants are in this cohort")
    site_path = fetch("cptac_coad_mutation_site.cbt", DOWNLOADS["cptac_coad_mutation_site.cbt"])
    site = pd.read_csv(site_path, sep="\t", index_col=0)
    n_patients = site.shape[1]
    kras = site.loc[[i for i in site.index if str(i).startswith("KRAS_")]]
    counts = (kras > 0).sum(axis=1).sort_values(ascending=False)
    counts = counts[counts > 0]
    out = pd.DataFrame({
        "variant": counts.index,
        "n_patients_mutated": counts.to_numpy(),
        "pct_of_cohort": (100 * counts / n_patients).to_numpy().round(2),
    })
    any_kras = int(((kras > 0).sum(axis=0) > 0).sum())
    log(f"  {any_kras}/{n_patients} patients ({100*any_kras/n_patients:.1f}%) have some "
        f"KRAS coding mutation, across {len(out)} different variants")
    log(out.to_string(index=False))
    STATS.update(kras_any=any_kras, kras_variants=len(out))
    return out


# Where I look for the independent-cohort spreadsheet. It's supplementary data
# from a journal, so it can't be downloaded automatically - you have to save it
# yourself. If it isn't there, this step is skipped and nothing else breaks.
FUDAN_XLSX_CANDIDATES = [
    os.path.join(RAW, "can-23-1878_supplementary_tables_1-5_suppst1-5.xlsx"),
    os.path.expanduser("~/Downloads/can-23-1878_supplementary_tables_1-5_suppst1-5.xlsx"),
]


def build_independent_proteomics() -> pd.DataFrame | None:
    """Tumour vs normal protein in a second, unrelated cohort, for checking.

    Source: supplementary Table S6B of "Comprehensive Proteogenomic Profiling
    Reveals the Molecular Characteristics of Colorectal Cancer at Distinct
    Stages of Progression", Cancer Research 2024;84:2888-2910 (Fudan cohort).
    That table has paired tumour and adjacent-normal proteomics for 28 patients,
    measured label-free rather than by TMT. Agreeing with CPTAC would be a much
    better sign than CPTAC agreeing with itself.
    """
    log("\n[6/6] CONTEXT: independent proteomics cohort (Fudan, Cancer Res 2024)")
    path = next((p for p in FUDAN_XLSX_CANDIDATES if os.path.exists(p)), None)
    if path is None:
        log("  spreadsheet not found, skipping. To include it, save "
            "can-23-1878_supplementary_tables_1-5_suppst1-5.xlsx into data/raw/")
        return None
    try:
        import openpyxl  # noqa: F401
    except ImportError:
        log("  openpyxl isn't installed, skipping this step "
            "(`conda install openpyxl` if you want it)")
        return None

    import re
    log(f"  reading {path}")
    wb = openpyxl.load_workbook(path, read_only=True)
    ws = wb["Supplementary Table S6B"]
    it = ws.iter_rows(values_only=True)
    next(it)                                   # title row
    hdr = [str(h).strip() if h is not None else "" for h in next(it)]
    ncol = max(i for i, h in enumerate(hdr) if h) + 1
    hdr = hdr[:ncol]

    data = {}
    for row in it:
        if not row or row[0] is None:
            continue
        vals = []
        for v in row[1:ncol]:
            try:
                vals.append(float(v))
            except (TypeError, ValueError):
                vals.append(np.nan)
        data[str(row[0]).strip().upper()] = vals
    df = pd.DataFrame.from_dict(data, orient="index", columns=hdr[1:ncol])
    log(f"  table S6B: {df.shape[0]} proteins x {df.shape[1]} samples")

    # Columns are named like L1_NAT and L1_Tumor, so I group them by the bit
    # before the underscore to find complete pairs.
    pref: dict[str, dict[str, str]] = {}
    for c in df.columns:
        m = re.match(r"^(.*)_(NAT|Tumor)$", c)
        if m:
            pref.setdefault(m.group(1), {})[m.group(2)] = c
    pairs = [(v["NAT"], v["Tumor"]) for v in pref.values() if "NAT" in v and "Tumor" in v]
    log(f"  {len(pairs)} patients have both a tumour and a normal sample")

    # This study writes missing values as 1e-05, so I treat that as not
    # detected rather than as a real measurement.
    FLOOR = 1e-05
    rows = []
    for g in df.index:
        ratios = []
        for n, t in pairs:
            a, b = df.at[g, n], df.at[g, t]
            if np.isfinite(a) and np.isfinite(b) and a > FLOOR and b > FLOOR:
                ratios.append(np.log2(b / a))
        if len(ratios) >= 5:
            pv = (stats.wilcoxon(ratios).pvalue
                  if len(ratios) >= 6 and any(x != 0 for x in ratios) else np.nan)
            rows.append((g, float(np.median(ratios)), len(ratios), pv))
    out = pd.DataFrame(rows, columns=["gene", "log2fc_fudan", "n_pairs_fudan", "pval_fudan"])
    log(f"  -> {len(out)} proteins with at least 5 usable pairs")
    STATS.update(n_fudan_pairs=len(pairs), n_fudan_genes=len(out))
    return out


# ---------------------------------------------------------------------- shared

def benjamini_hochberg(p: np.ndarray) -> np.ndarray:
    """BH-adjusted p-values (FDR), skipping NaNs."""
    q = np.full(len(p), np.nan)
    ok = ~np.isnan(p)
    vals = p[ok]
    n = len(vals)
    if n == 0:
        return q
    order = np.argsort(vals)
    ranked = vals[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.clip(ranked, 0, 1)
    q[ok] = adj
    return q


PROVENANCE = """# Where the colon data came from

Everything in this folder was made by `prepare_coad_data.py` in the repo root,
from real downloads. Nothing here is made up or simulated. Running that script
again rebuilds all of it.

## data/coad_transcriptomics.tsv
* **Study** TCGA Colon Adenocarcinoma (TCGA-COAD), TCGA Network, *Nature*
  2012;487:330-337.
* **Files** `TCGA-COAD.star_tpm.tsv.gz` and
  `gencode.v36.annotation.gtf.gene.probemap`, from the UCSC Xena GDC hub at
  `https://gdc-hub.s3.us-east-1.amazonaws.com/download/`.
* **What the source contains** log2(TPM+1) values, STAR quantification, GENCODE
  v36, rows labelled with Ensembl gene IDs.
* **Samples** {n_tumor} primary tumours (barcode type `-01`) and {n_normal}
  adjacent normals (`-11`). I left out the one recurrent (`-02`) and one
  metastatic (`-06`) sample.
* **What I did** mapped Ensembl IDs to gene symbols, keeping the
  highest-expressed ID where a symbol had several; dropped genes with mean
  log2(TPM+1) below 1 in both groups; took `log2fc` = mean(tumour) -
  mean(normal), which is a real log2 fold change because the values are already
  logged; and used a two-sided Mann-Whitney U test for `pval`, since these are
  different patients and so not paired.
* **Rows** {n_tx}

## data/coad_proteomics.tsv
* **Study** CPTAC-2 Prospective Colon Cancer, Vasaikar et al., *Cell*
  2019;177:1035-1049.e19 (PMID 31031003, doi 10.1016/j.cell.2019.03.030).
* **Files** `...PNNL_Tumor_TMT_UnsharedLogRatio.cct` and
  `...PNNL_Normal_TMT_UnsharedLogRatio.cct`, from the LinkedOmics CPTAC-COAD
  cohort.
* **What the source contains** gene-level TMT log2 ratios against the reference
  channel, using unshared peptides only. Rows are gene symbols.
* **Samples** {n_pairs} patients with both a tumour and an adjacent-normal
  channel, so this comparison is properly paired.
* **What I did** took tumour minus normal within each patient, used the mean of
  those differences as `log2fc`, and a Wilcoxon signed-rank test on them for
  `pval`. Genes with fewer than {min_pairs} usable pairs were dropped.
* **Check** my log2FC matches the study's own `Tumor_Normal_log2FC.cct` at
  Pearson r = {prot_r:.4f} across {n_shared} shared genes. That's how I know
  this layer is right.
* **Rows** {n_pr}

## data/coad_genomics.tsv
* **Study** the whole-exome data from the same CPTAC-2 colon study, so the same
  patients as the protein layer. WUSM GATK pipeline, gene-level calls.
* **File** `...WUSM__Mutation__GAIIx__03_01_2017__BCM__Gene__GATK_Pipeline.cbt`
  from the LinkedOmics CPTAC-COAD cohort.
* **What the source contains** a 1 if that patient's tumour has a somatic
  mutation in that gene, 0 if not. {n_patients} patients.
* **Gene lengths** longest CDS per gene symbol, from UCSC hg38 `refFlat.txt.gz`.
* **What I did** `mut_freq` is the fraction of patients mutated. Then I worked
  out a background rate of mu = {mu:.3e} mutated-gene events per coding base
  per patient across all genes with a known CDS length, converted it to a
  per-gene chance of being hit as p = 1 - exp(-mu x CDS length), and used a
  one-sided binomial test on the observed number of mutated patients.
  `neglog10p` is -log10 of that p-value.
* **Why bother with the length correction** without it the top of the list is
  TTN (56.6% of patients), MUC16 (36.8%) and OBSCN (27.4%), all of which
  outrank or match KRAS (33.0%) just by being long. This is a length
  correction in the same spirit as MutSigCV, not a proper driver-gene method,
  and I don't present it as one.
* **Rows** {n_gw}

## Quick check that the mutation file is what it says it is
Straight out of the CPTAC file, the gene-level frequencies are APC 77.4%,
TP53 52.8%, KRAS 33.0%, BRAF 17.0%. That's the standard colon adenocarcinoma
picture, which is reassuring.

## data/coad_cptac_matched_rna_protein_corr.tsv (context, not used for scoring)
Per-gene Spearman correlation between CPTAC tumour RNA
(`...UNC__RNAseq...BCM_RSEM_UpperQuartile_log2.cct`) and CPTAC tumour protein,
across the {n_matched_samples} tumours measured on both. {n_corr_genes} genes
with at least 20 usable samples, median Spearman r = {median_corr:.3f}. This is
only here so I can say what the cross-cohort design costs, using numbers from
inside one cohort.

## data/coad_kras_variant_spectrum.tsv (context, not used for scoring)
From the site-level calls in
`...WUSM__Mutation__GAIIx__03_01_2017__BCM__Site__GATK_Pipeline.cbt`.
{kras_any} of {n_patients} patients have a KRAS coding mutation, across
{kras_variants} different variants.

## data/coad_independent_proteomics_fudan.tsv (context, not used for scoring)
Supplementary Table S6B of "Comprehensive Proteogenomic Profiling Reveals the
Molecular Characteristics of Colorectal Cancer at Distinct Stages of
Progression", *Cancer Research* 2024;84:2888-2910. A separate Chinese (Fudan)
cohort with paired tumour and adjacent-normal proteomics, measured label-free
instead of by TMT. {n_fudan_pairs} complete pairs, {n_fudan_genes} proteins
with at least 5 usable pairs. Missing values are written as 1e-05 in that
table and I treat them as not detected. I use it to check the CPTAC protein
layer against a completely separate group of patients.

## About the two different cohorts
The RNA layer is TCGA-COAD; the protein and mutation layers are CPTAC-2 colon.
**No patient is in both.** So everything is joined at the gene level and I
can't say anything about individual patients. CPTAC does have RNA-seq on its
own 106 patients, but only for the tumours, with no adjacent-normal RNA, so it
can't give a tumour vs normal RNA comparison. That's the reason the RNA layer
had to come from TCGA.
"""


def main() -> None:
    os.makedirs(RAW, exist_ok=True)
    os.makedirs(OUT, exist_ok=True)

    tx = build_transcriptomics()
    pr = build_proteomics()
    gw = build_genomics()

    tx_cols = ["gene", "log2fc", "pval", "qval", "mean_tumor_log2tpm",
               "mean_normal_log2tpm", "n_tumor", "n_normal"]
    pr_cols = ["gene", "log2fc", "pval", "qval", "n_pairs"]
    gw_cols = ["gene", "neglog10p", "pval", "qval", "mut_freq", "n_mutated",
               "cds_len", "n_patients"]

    tx[tx_cols].to_csv(f"{OUT}/coad_transcriptomics.tsv", sep="\t", index=False)
    pr[pr_cols].to_csv(f"{OUT}/coad_proteomics.tsv", sep="\t", index=False)
    gw[gw_cols].to_csv(f"{OUT}/coad_genomics.tsv", sep="\t", index=False)

    log("\nWrote:")
    for f in ("coad_transcriptomics.tsv", "coad_proteomics.tsv", "coad_genomics.tsv"):
        log(f"  {OUT}/{f}  ({os.path.getsize(f'{OUT}/{f}')/1e3:.0f} kB)")

    corr = build_matched_rna_protein_correlation()
    corr.to_csv(f"{OUT}/coad_cptac_matched_rna_protein_corr.tsv", sep="\t", index=False)
    log(f"  {OUT}/coad_cptac_matched_rna_protein_corr.tsv")

    kras = build_kras_variant_spectrum()
    kras.to_csv(f"{OUT}/coad_kras_variant_spectrum.tsv", sep="\t", index=False)
    log(f"  {OUT}/coad_kras_variant_spectrum.tsv")

    fudan = build_independent_proteomics()
    if fudan is not None:
        fudan.to_csv(f"{OUT}/coad_independent_proteomics_fudan.tsv", sep="\t", index=False)
        log(f"  {OUT}/coad_independent_proteomics_fudan.tsv")
    else:
        STATS.setdefault("n_fudan_pairs", "n/a")
        STATS.setdefault("n_fudan_genes", "n/a")

    overlap = set(tx["gene"]) & set(pr["gene"]) & set(gw["gene"])
    log(f"\nGenes present in all three layers: {len(overlap)}")

    with open(f"{OUT}/PROVENANCE.md", "w") as fh:
        fh.write(PROVENANCE.format(**STATS))
    log(f"Wrote {OUT}/PROVENANCE.md")


if __name__ == "__main__":
    main()
