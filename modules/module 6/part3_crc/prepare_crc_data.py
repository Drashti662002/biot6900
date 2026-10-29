"""
prepare_crc_data.py

This builds the four files that the Part 3 notebook reads, for colorectal
cancer instead of Alzheimer's.

Why I needed it: the Canvas data.zip is built around an Alzheimer's score
(PGS000334), and its genotypes.csv only holds those 22 SNPs, so a different
disease's score has nothing to match against. My disease since Module 2 has been
colon adenocarcinoma, so I rebuilt the same four files myself from the same
sources the professor used: a PGS Catalog score and 1000 Genomes phase 3. I
also rebuilt the simulated phenotype with the professor's recipe (from the
Canvas data/README.md), changing nothing except the score.

Run it once from this folder (needs bcftools on PATH and an internet connection):

    python prepare_crc_data.py

It writes:

    data/genotypes.csv   sample_id + one column per score SNP, ALT dosage 0/1/2
    data/variants.csv    variant, ref, alt (from the 1000 Genomes VCF)
    data/score.csv       variant, effect_allele, other_allele, weight (as published)
    data/meta.csv        sample_id, ancestry, age, phenotype (phenotype SIMULATED)
    data/simulation_log.csv   one row per tag: the causal SNP it got, and its r2 per ancestry
    data/README.md       sources, cohort, simulation, checksums

Sources (all open access):

1. Score - PGS Catalog PGS000148 (CRC63; Jeon J et al., Gastroenterology
   2018;154:2152-2164, doi:10.1053/j.gastro.2018.02.021). 63 SNPs, weights are
   per-allele log(OR) estimated in European-ancestry samples only. I use the
   harmonized GRCh37 file.

2. Genotypes - 1000 Genomes phase 3 (release 20130502, GRCh37), read
   region-by-region from the EBI FTP with bcftools, so nothing chromosome-sized
   is downloaded.

3. People - the same 2,336 samples, ancestry labels and (simulated) ages as the
   Canvas meta.csv, so my cohort is identical to the class cohort and only the
   score and the phenotype differ.
"""

from __future__ import annotations

import gzip
import hashlib
import os
import subprocess
import sys
import urllib.request

import numpy as np
import pandas as pd

# ------------------------------------------------------------------- settings

RAW = "data/raw"          # downloads and cached regions, git ignores this folder
OUT = "data"              # the four small files the notebook reads
CLASS_META = "../data/meta.csv"   # the Canvas cohort (sample ids, ancestry, age)

PGS_ID = "PGS000148"
SCORE_URL = ("https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/"
             f"{PGS_ID}/ScoringFiles/Harmonized/{PGS_ID}_hmPOS_GRCh37.txt.gz")
KG = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
      "ALL.chr{c}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")

# The professor's simulation recipe (Canvas data/README.md), unchanged
SEED = 6900
WINDOW = 250_000          # causal searched within +/- 250 kb of the tag
R2_MIN = 0.6              # r2 with the tag in EUR
EUR_MAF_MIN = 0.05        # causal must be common in EUR
AGE_COEF = 0.35           # liability += 0.35 * (age - 60) / 8
NOISE_SD = 2.09
CASE_FRACTION = 0.30      # top 30% of liability across the whole cohort are cases

POPS = ["EUR", "AFR", "EAS", "SAS", "AMR"]
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def log(msg: str) -> None:
    print(msg, flush=True)


def sha8(path: str) -> str:
    return hashlib.sha256(open(path, "rb").read()).hexdigest()[:8]


# ---------------------------------------------------------------------- score

def load_score() -> pd.DataFrame:
    path = os.path.join(RAW, f"{PGS_ID}_hmPOS_GRCh37.txt.gz")
    if not os.path.exists(path):
        urllib.request.urlretrieve(SCORE_URL, path)
    with gzip.open(path, "rt") as fh:
        s = pd.read_csv(fh, sep="\t", comment="#", dtype=str)
    s = s.rename(columns={"hm_rsID": "variant", "effect_weight": "weight"})
    s["chr"] = s["hm_chr"]
    s["pos"] = s["hm_pos"].astype(int)
    s["weight"] = s["weight"].astype(float)
    s["effect_allele"] = s["effect_allele"].fillna("")
    s["other_allele"] = s["other_allele"].fillna("")
    return s[["variant", "chr", "pos", "effect_allele", "other_allele", "weight"]]


# ------------------------------------------------------------------ genotypes

def windows(score: pd.DataFrame) -> list[tuple[str, int, int]]:
    """+/- WINDOW around every tag, merged where they overlap on a chromosome."""
    out = []
    for c, g in score.sort_values("pos").groupby("chr", sort=False):
        lo = hi = None
        for p in g["pos"]:
            a, b = max(1, p - WINDOW), p + WINDOW
            if lo is None or a > hi:
                if lo is not None:
                    out.append((c, lo, hi))
                lo, hi = a, b
            else:
                hi = max(hi, b)
        out.append((c, lo, hi))
    return out


def fetch_region(c: str, lo: int, hi: int, tag_pos: list[int], samples: str) -> str:
    """Biallelic SNPs with MAF >= 1% in the cohort, plus the tag positions at any MAF."""
    cache = os.path.join(RAW, "regions", f"chr{c}_{lo}_{hi}.npz")
    if os.path.exists(cache):
        return cache
    keep_tags = " || ".join(f"POS=={p}" for p in tag_pos)
    view = ["bcftools", "view", "-S", samples, "-r", f"{c}:{lo}-{hi}",
            "-m2", "-M2", "-v", "snps", "-Ou", KG.format(c=c)]
    filt = ["bcftools", "view", "-i", f"MAF>=0.01 || {keep_tags}", "-Ou"]
    query = ["bcftools", "query", "-f", "%POS\t%ID\t%REF\t%ALT[\t%GT]\n"]
    order = subprocess.run(["bcftools", "view", "-h", "-S", samples, KG.format(c=c)],
                           capture_output=True, text=True, check=True, cwd=RAW)
    order = order.stdout.strip().splitlines()[-1].split("\t")[9:]
    # the EBI connection sometimes hangs, so give each try 5 minutes and retry
    pipeline = " | ".join(" ".join(f"'{a}'" for a in cmd) for cmd in (view, filt, query))
    for attempt in range(4):
        try:
            res = subprocess.run(pipeline, shell=True, capture_output=True, text=True,
                                 cwd=RAW, timeout=300, check=True)
            break
        except (subprocess.TimeoutExpired, subprocess.CalledProcessError):
            log(f"    retry {attempt + 1} for chr{c}:{lo}-{hi}")
    else:
        raise RuntimeError(f"bcftools failed on chr{c}:{lo}-{hi}")
    pos, ids, ref, alt, dos = [], [], [], [], []
    for line in res.stdout.splitlines():
        f = line.split("\t")
        pos.append(int(f[0])); ids.append(f[1]); ref.append(f[2]); alt.append(f[3])
        dos.append(np.fromiter((g.count("1") for g in f[4:]), dtype=np.int8, count=len(f) - 4))
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    np.savez_compressed(cache, pos=np.array(pos), ids=np.array(ids), ref=np.array(ref),
                        alt=np.array(alt), G=np.vstack(dos), samples=np.array(order))
    return cache


# ----------------------------------------------------------------- phenotype

def corr_rows(x: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """Pearson r between one dosage vector x and every row of Y."""
    x = x - x.mean()
    Yc = Y - Y.mean(axis=1, keepdims=True)
    den = np.sqrt((x ** 2).sum() * (Yc ** 2).sum(axis=1))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, Yc @ x / den, 0.0)


def main() -> None:
    os.makedirs(RAW, exist_ok=True)
    meta = pd.read_csv(CLASS_META)[["sample_id", "ancestry", "age"]]
    samples = os.path.abspath(os.path.join(RAW, "samples_2336.txt"))
    meta["sample_id"].to_csv(samples, index=False, header=False)

    score = load_score()
    log(f"{PGS_ID}: {len(score)} variants")

    # 1. pull every region, find each tag by position and alleles
    tag_rows, regions = {}, {}
    for c, lo, hi in windows(score):
        in_win = score[(score.chr == c) & score.pos.between(lo, hi)]
        log(f"  chr{c}:{lo:,}-{hi:,}  ({len(in_win)} tags)")
        z = np.load(fetch_region(c, lo, hi, list(in_win.pos), samples))
        idx = pd.Index(z["samples"]).get_indexer(meta["sample_id"])
        assert (idx >= 0).all(), "a cohort sample is missing from the VCF"
        R = dict(pos=z["pos"], ids=z["ids"], ref=z["ref"], alt=z["alt"], G=z["G"][:, idx])
        regions[(c, lo, hi)] = R
        for _, t in in_win.iterrows():
            hits = np.where(R["pos"] == t.pos)[0]
            ok = [h for h in hits
                  if (t.effect_allele == "" or t.effect_allele in (R["ref"][h], R["alt"][h]))
                  and (t.other_allele == "" or t.other_allele in (R["ref"][h], R["alt"][h]))]
            tag_rows[t.variant] = ((c, lo, hi), ok[0]) if ok else None

    missing = [v for v, r in tag_rows.items() if r is None]
    log(f"tags found in 1000 Genomes: {len(score) - len(missing)} of {len(score)}"
        + (f"  missing: {missing}" if missing else ""))
    found = [v for v in score.variant if tag_rows.get(v) is not None]

    # 2. the four notebook files
    geno = pd.DataFrame({v: regions[tag_rows[v][0]]["G"][tag_rows[v][1]] for v in found})
    geno.insert(0, "sample_id", meta["sample_id"].values)
    variants = pd.DataFrame({"variant": found,
                             "ref": [regions[tag_rows[v][0]]["ref"][tag_rows[v][1]] for v in found],
                             "alt": [regions[tag_rows[v][0]]["alt"][tag_rows[v][1]] for v in found]})
    sc = score.set_index("variant").loc[found].reset_index()
    sc = sc[["variant", "effect_allele", "other_allele", "weight"]]
    # one published row has no other_allele; fill it from the VCF like the class file does
    # (a row with no effect_allele is left blank on purpose: the notebook should drop it)
    blank = (sc["other_allele"] == "") & (sc["effect_allele"] != "")
    v_idx = variants.set_index("variant")
    sc.loc[blank, "other_allele"] = [
        v_idx.at[v, "alt"] if e == v_idx.at[v, "ref"] else v_idx.at[v, "ref"]
        for v, e in zip(sc.loc[blank, "variant"], sc.loc[blank, "effect_allele"])]
    filled_other = list(sc.loc[blank, "variant"])
    no_effect = list(sc.loc[sc["effect_allele"] == "", "variant"])

    # 3. simulate the phenotype, the professor's recipe
    rng = np.random.default_rng(SEED)
    anc = meta["ancestry"].to_numpy()
    eur = anc == "EUR"
    genetic = np.zeros(len(meta))
    sim = []
    for _, t in sc.iterrows():
        key, h = tag_rows[t.variant]
        R = regions[key]
        G = R["G"].astype(float)
        tag = G[h].copy()
        # no published effect allele: the disease still has this locus, so for the
        # simulation only, take the effect allele as the one that is not other_allele
        eff_sim = t.effect_allele or (R["alt"][h] if t.other_allele == R["ref"][h] else R["ref"][h])
        if eff_sim == R["ref"][h]:
            tag = 2 - tag                                    # count the effect allele
        r_eur = corr_rows(tag[eur], G[:, eur])
        f_eur = G[:, eur].mean(axis=1) / 2
        maf_eur = np.minimum(f_eur, 1 - f_eur)
        near = np.abs(R["pos"] - R["pos"][h]) <= WINDOW
        cand = np.where(near & (r_eur ** 2 > R2_MIN) & (maf_eur >= EUR_MAF_MIN)
                        & (np.arange(len(r_eur)) != h))[0]
        k = rng.choice(cand) if len(cand) else h             # no neighbour -> own causal
        sign = np.sign(r_eur[k]) if k != h else 1.0
        causal = G[k] if k != h else tag
        genetic += t.weight * sign * causal
        r2 = {p: float(np.corrcoef(tag[anc == p], G[k, anc == p])[0, 1] ** 2) if k != h
              else 1.0 for p in POPS}
        sim.append({"tag": t.variant, "weight": t.weight,
                    "causal": R["ids"][k] if R["ids"][k] != "." else f"{key[0]}:{R['pos'][k]}",
                    "causal_pos": int(R["pos"][k]), "own_causal": k == h,
                    "n_candidates": len(cand), **{f"r2_{p}": round(r2[p], 3) for p in POPS}})
    genetic = (genetic - genetic.mean()) / genetic.std()
    age = meta["age"].to_numpy()
    liability = genetic + AGE_COEF * (age - 60) / 8 + rng.normal(0, NOISE_SD, len(meta))
    phenotype = (liability >= np.quantile(liability, 1 - CASE_FRACTION)).astype(int)
    meta_out = meta.assign(phenotype=phenotype)

    # 4. write
    paths = {"genotypes.csv": geno, "variants.csv": variants, "score.csv": sc, "meta.csv": meta_out}
    for name, df in paths.items():
        df.to_csv(os.path.join(OUT, name), index=False)
    sim = pd.DataFrame(sim)
    sim.to_csv(os.path.join(OUT, "simulation_log.csv"), index=False)

    # quick checks against the notebook's own harmonization rule
    eff = sc.set_index("variant")["effect_allele"]
    keep = int((eff == v_idx["alt"].loc[eff.index]).sum())
    flip = int((eff == v_idx["ref"].loc[eff.index]).sum())
    ambiguous = [v for v in found if {v_idx.at[v, "ref"], v_idx.at[v, "alt"]}
                 in ({"A", "T"}, {"C", "G"})]
    mean_r2 = {p: round(sim[f"r2_{p}"].mean(), 2) for p in POPS}
    log(f"harmonization: keep {keep}, flip {flip}, drop {len(found) - keep - flip}")
    log(f"strand-ambiguous (A/T or C/G): {ambiguous}")
    log(f"own causal (no neighbour): {list(sim.loc[sim.own_causal, 'tag'])}")
    log(f"mean tag-causal r2: {mean_r2}")
    log(f"cases: {phenotype.mean():.3f} overall; "
        + ", ".join(f"{p} {phenotype[anc == p].mean():.3f}" for p in POPS))
    write_readme(score, found, missing, filled_other, no_effect, keep, flip, ambiguous, sim, mean_r2,
                 meta_out, {f: sha8(os.path.join(OUT, f)) for f in paths})


def write_readme(score, found, missing, filled_other, no_effect, keep, flip, ambiguous, sim, mean_r2,
                 meta, hashes) -> None:
    counts = meta["ancestry"].value_counts()
    case = (100 * meta.groupby("ancestry")["phenotype"].mean()).round(0).astype(int)
    own = list(sim.loc[sim.own_causal, "tag"])
    lines = f"""# BIOT 6900 · Module 6, Part 3 cohort: colorectal cancer

> **THE PHENOTYPE IS SIMULATED.** 1000 Genomes has no disease labels. The genotypes, alleles and score weights are real; `phenotype` in `meta.csv` was generated by the simulation below. Age is also simulated (it is the Canvas cohort's age column). AUC and calibration results illustrate a real phenomenon; they are not measurements of real colorectal cancer risk.

Built by `../prepare_crc_data.py`. Same four files and columns as the Canvas `data.zip`, so the Part 3 notebook reads them unchanged.

## Sources

- **Score:** PGS Catalog **{PGS_ID}** (CRC63; Jeon J et al., *Gastroenterology* 2018;154:2152-2164, doi:10.1053/j.gastro.2018.02.021). {len(score)} SNPs, weights are per-allele log(OR) from a European-ancestry-only GWAS. Harmonized scoring file, **GRCh37** positions. Published EUR validation AUROC 0.59 (PGS Catalog PPM000450, with age, family history, study and endoscopy history as covariates).
  - One row ({", ".join(filled_other) or "none"}) publishes no `other_allele`; I filled it with whichever of the 1000 Genomes REF/ALT is not the effect allele. Effect alleles and weights are unchanged.
  - One row ({", ".join(no_effect) or "none"}) publishes **no `effect_allele`** (only `other_allele`). I left it blank, so the notebook's keep/flip/drop rule drops it. The simulated disease still includes this locus (taking the effect allele as the non-`other_allele` one), because the biology is real even when the scoring file is broken.
- **Genotypes:** **1000 Genomes phase 3** (release 20130502, **GRCh37**), biallelic SNPs, read region by region from the EBI FTP with bcftools. Tags matched by chromosome, position and alleles.
- **People:** the same 2,336 samples, super-population labels and ages as the Canvas `meta.csv`.

## Cohort

| Ancestry | N | % cases (simulated) |
|---|---|---|
""" + "\n".join(f"| {p} | {counts[p]} | {case[p]} |" for p in POPS) + f"""
| **Total** | **{len(meta)}** | {round(100 * meta.phenotype.mean())} |

Score coverage: **{len(found)} of {len(score)}** variants present{(" (missing: " + ", ".join(missing) + ", which has no record at its GRCh37 position in 1000 Genomes phase 3)") if missing else ""}. Harmonization gives keep {keep}, flip {flip}, drop {len(found) - keep - flip}. Strand-ambiguous SNPs (A/T or C/G), kept because the PGS Catalog file is on the + strand like the VCF: {", ".join(ambiguous) or "none"}.

## How the phenotype was simulated

The professor's recipe from the Canvas README, with the same numbers. Random seed **{SEED}**.

1. For each score (tag) SNP, a hidden **causal** SNP was chosen at random among 1000 Genomes SNPs within ±250 kb with **r² > {R2_MIN} with the tag in EUR** and EUR minor-allele frequency ≥ 5%. Tags with no such neighbour are their own causal: {", ".join(own) or "none"}.
2. Each causal SNP got a true per-allele effect equal to its tag's published weight, signed by the direction of its EUR LD with the tag.
3. Liability = genetic score from the causal SNPs (scaled to SD 1) + {AGE_COEF} × (age − 60)/8 + Normal(0, {NOISE_SD}²) noise.
4. The top {int(CASE_FRACTION * 100)}% of liability across the whole cohort are cases. **The same rule applies to every ancestry.**

I did **not** retune the noise SD for this score; it is the professor's value. Mean r² between tags and their causal SNPs: """ + ", ".join(f"{p} {mean_r2[p]}" for p in POPS) + f""". Per-tag details in `simulation_log.csv`.

## Limitations

- Causal SNPs are picked from SNPs in LD with each tag **in Europeans**, so by construction the score tags its causal SNPs best in EUR. That is how the portability drop arises here, from real 1000 Genomes LD, and it is also why the drop's size depends on this design.
- PGS000148 weights come from a joint model, so some loci carry several SNPs with opposite-sign weights (e.g. the FADS1/2 SNPs on chr11). Each tag still gets its own causal SNP, so in those loci the simulated causal architecture is more complicated than the real one probably is.

## Checksums (first 8 hex characters of SHA-256, as printed by the notebook)

""" + " · ".join(f"{f} `{h}`" for f, h in hashes.items()) + "\n"
    open(os.path.join(OUT, "README.md"), "w").write(lines)


if __name__ == "__main__":
    sys.exit(main())
