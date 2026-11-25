# Lab 6 (Module 6 Part 3, CRC PRS) evidence inventory

Base path: `/Users/drashtibhanushali/biot6900/modules/module 6/`
- NB = `part3_crc/BIOT6900_Module6_Part3_Graded_CRC.ipynb` (cell indices 0-based; code cells execution_count 1-10 sequential, i.e. Restart & Run All; kernel "Python 3 (biot6900)")
- REPORT = `part3_crc/module 6_report.docx` (pandoc plain text)
- DREADME = `part3_crc/data/README.md`; TREADME = `README.md` (module root)
- Plots: NB cell 14 -> PRS density histogram by ancestry (saved `part3_crc/figures/prs_by_ancestry.png`; report Figure 1). NB cell 19 -> calibration curve by ancestry (saved `part3_crc/figures/calibration_by_ancestry.png`; report Figure 2). No other images in NB.

## 1. PGS000148 harmonisation

| Item | Value | Provenance |
|---|---|---|
| Published SNPs | 63 | raw scoring file header `#variants_number=63` (`part3_crc/data/raw/PGS000148_hmPOS_GRCh37.txt.gz`); DREADME; REPORT §3.1 |
| In score.csv / genotype panel | 62 | NB cell 8 output "score size : 62 variants" |
| Used | 61 | NB cell 10 output "used 61  flipped 28  dropped 1" |
| Flipped (effect allele = REF) | 28 | NB cell 10 output; list: rs10911251, rs72647484, rs10936599, rs2450115, rs10505477, rs6983267, rs7014346, rs719725, rs10904849, rs10795668, rs4919687, rs1535, rs174537, rs174550, rs4246215, rs3802842, rs10774214, rs11064437, rs11169552, rs7136702, rs1957636, rs4779584, rs79900961, rs9929218, rs4939827, rs10411210, rs6066825, rs6061231 |
| Kept as-is | 33 | DREADME "keep 33, flip 28, drop 1" |
| Dropped | rs7229639 | NB cell 10 output "dropped: ['rs7229639']" |
| Why dropped | Published row has blank effect_allele (only other_allele = G, weight 0.023); blank string matches neither REF (A) nor ALT (G), so keep/flip/drop drops it | raw file row `rs7229639 18 <blank> G 0.023 0.09 ... 18 46450976`; `data/score.csv` row `rs7229639,,G,0.023`; `data/variants.csv` `rs7229639,A,G`; REPORT §3.1 ("The effect_allele in the published row is blank..."); TREADME "What didn't work" |
| rs7229639 in simulation | Kept in simulated disease (effect allele taken as non-other_allele) | DREADME Sources; `prepare_crc_data.py` (eff_sim); `simulation_log.csv` row rs7229639 (causal 18:46452637, r2 EUR 0.944, AFR 0.23, EAS 0.03, SAS 0.942, AMR 0.24) |
| Not found | rs202110856 (chr5:96129873 GRCh37, G/C, weight 0.02, published effect-allele freq 0.99) | raw scoring file row; DREADME "no record at its GRCh37 position in 1000 Genomes phase 3"; REPORT §3.1 |
| Strand-ambiguous kept | rs812481, rs35360328, rs16941835 | DREADME; REPORT §3.1 |
| Own-causal tags (no LD neighbour) | rs704017, rs3217810, rs11064437, rs79900961, rs12603526 | DREADME step 1; `simulation_log.csv` own_causal=True |

## 2. Per-ancestry table

Sources: n / mean PRS / SD / % above cutoff / % cases = NB cell 14 text output; AUC = NB cell 16 output; O/E, predicted, observed = NB cell 19 output; mean tag-causal r² = DREADME (not printed in NB); recomputed column means of `data/simulation_log.csv` r2_* columns shown in brackets (AFR and EAS each have 1 NaN, skipped).

Cutoff: "pooled top-10% cutoff: PRS > 1.282" = `np.quantile(meta["prs"], 0.90)` of the cohort-standardised (z-scored, pooled mean 0/SD 1) PRS (NB cells 12, 14). 1.282 is the empirical 90th percentile of the pooled z-score (coincides numerically with the normal 90th pct).

| Ancestry | n | cases | mean PRS (z) | SD PRS | % > 1.282 | % cases (sim) | AUC | mean tag-causal r² | mean predicted | observed | O/E | pred range |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| EUR | 500 | 150 | 0.16 | 0.98 | 12.2 | 30.0 | 0.685 | 0.83 [0.826] | 0.300 | 0.300 | 1.00 | 0.04-0.79 |
| AFR | 500 | 133 | -0.09 | 1.00 | 8.4 | 26.6 | 0.567 | 0.45 [0.452] | 0.268 | 0.266 | 0.99 | 0.03-0.68 |
| EAS | 500 | 141 | -0.17 | 0.96 | 7.2 | 28.2 | 0.689 | 0.76 [0.757] | 0.256 | 0.282 | 1.10 | 0.05-0.73 |
| SAS | 489 | 186 | 0.30 | 0.93 | 14.9 | 38.0 | 0.606 | 0.77 [0.772] | 0.318 | 0.380 | 1.20 | 0.07-0.81 |
| AMR | 347 | 91 | -0.28 | 1.02 | 6.3 | 26.2 | 0.670 | 0.74 [0.737] | 0.244 | 0.262 | 1.08 | 0.04-0.72 |
| Overall | 2336 | 701* | 0 (by construction) | 1 | 10.0 (report only) | 30.0 | 0.647 | - | - | - | - | - |

*701 = sum of per-ancestry cases from NB cell 16; `meta.csv` overall phenotype mean 0.3001. Overall 10.0% above cutoff appears in REPORT table only (true by definition, not printed by NB).

## 3. EUR-fitted logistic model

NB cell 19 output: "EUR-fitted model: intercept -1.054, slope per SD 0.720". sklearn `LogisticRegression()` defaults (L2, C=1) fit on EUR z-PRS only -> applied to all. REPORT §3 Discrimination vs calibration quotes the same.

## 4. Per-ancestry quintile calibration — NUMBERS ARE SAVED AS TEXT

Verdict: the actual quintile numbers are printed in NB cell 19 stream output (not only the plot). Quintiles are within-ancestry quintiles of EUR-fitted predicted risk (`calibration_points(y, p, bins=5)`, NB cell 18), mean predicted vs observed case fraction per bin, 2 d.p. They are also trivially recomputable from saved data + saved code (deterministic: z-score PRS, sklearn logistic fit, quantile bins) without new modelling — not run.

| Ancestry | | Q1 | Q2 | Q3 | Q4 | Q5 |
|---|---|---|---|---|---|---|
| EUR | pred | 0.13 | 0.22 | 0.29 | 0.36 | 0.51 |
| EUR | obs | 0.13 | 0.23 | 0.26 | 0.36 | 0.52 |
| AFR | pred | 0.11 | 0.18 | 0.25 | 0.32 | 0.47 |
| AFR | obs | 0.16 | 0.26 | 0.31 | 0.26 | 0.34 |
| EAS | pred | 0.11 | 0.17 | 0.23 | 0.31 | 0.45 |
| EAS | obs | 0.10 | 0.24 | 0.23 | 0.34 | 0.50 |
| SAS | pred | 0.15 | 0.23 | 0.30 | 0.38 | 0.53 |
| SAS | obs | 0.32 | 0.32 | 0.30 | 0.39 | 0.58 |
| AMR | pred | 0.10 | 0.16 | 0.22 | 0.29 | 0.45 |
| AMR | obs | 0.11 | 0.20 | 0.25 | 0.29 | 0.46 |

Expected features check: AFR curve flat (obs 0.16-0.34 vs pred 0.11-0.47; non-monotone, Q4 0.26 < Q3 0.31); AFR top quintile over-predicted 0.47 vs 0.34; SAS Q1 obs/pred 0.32/0.15 = 2.1x, Q2 0.32/0.23 = 1.4x; SAS O/E 1.20; AFR O/E 0.99; AFR AUC 0.567. Bin sizes ~70 (AMR) to ~100 (n=500 groups) per REPORT; not printed in NB.

## 5. Age vs age+PRS

NB cell 21 output: "AUC age only: 0.579", "AUC age + PRS: 0.666" (5-fold `cross_val_score`, pooled cohort, sklearn default StratifiedKFold on case status, no shuffle, not stratified by ancestry, LogisticRegression, z-scored age and PRS). REPORT §3.6 and Q4: gain 0.087. Comparison quoted in REPORT/DREADME: PGS Catalog PPM000450 EUR AUROC 0.59 (with age, family history, study, endoscopy history covariates) — external figure, not computed here.

## 6. Simulated vs real

Real:
- Genotypes: 1000 Genomes phase 3 (release 20130502, GRCh37), biallelic SNPs, bcftools region reads (DREADME Sources; `prepare_crc_data.py` docstring).
- Weights/alleles: PGS Catalog PGS000148 (CRC63; Jeon J et al., Gastroenterology 2018;154:2152-2164, doi:10.1053/j.gastro.2018.02.021), harmonised GRCh37 file, per-allele log(OR) from European-ancestry-only GWAS (DREADME).
- People: same 2,336 samples and super-population labels as the Canvas meta.csv (DREADME).
- LD between tags and causal SNPs: real 1000 Genomes LD (DREADME Limitations).

Simulated (quotes):
- DREADME banner: "THE PHENOTYPE IS SIMULATED. 1000 Genomes has no disease labels. The genotypes, alleles and score weights are real; `phenotype` in `meta.csv` was generated by the simulation below. Age is also simulated (it is the Canvas cohort's age column). AUC and calibration results illustrate a real phenomenon; they are not measurements of real colorectal cancer risk."
- DREADME recipe (seed 6900): "1. For each score (tag) SNP, a hidden causal SNP was chosen at random among 1000 Genomes SNPs within ±250 kb with r² > 0.6 with the tag in EUR and EUR minor-allele frequency ≥ 5%. ... 2. Each causal SNP got a true per-allele effect equal to its tag's published weight, signed by the direction of its EUR LD with the tag. 3. Liability = genetic score from the causal SNPs (scaled to SD 1) + 0.35 × (age − 60)/8 + Normal(0, 2.09²) noise. 4. The top 30% of liability across the whole cohort are cases. The same rule applies to every ancestry."
- `prepare_crc_data.py` constants: SEED=6900, WINDOW=250_000, R2_MIN=0.6, EUR_MAF_MIN=0.05, AGE_COEF=0.35, NOISE_SD=2.09, CASE_FRACTION=0.30.
- NB cell 2 (student note): "The phenotype and age are simulated, as in the Canvas data. AUC and calibration results are illustrative."
- REPORT Setup: "The phenotype is simulated, and so is the age effect. I followed the recipe in the Canvas README and didn't change any of the settings (seed 6900, ±250 kb window, r² > 0.6 in EUR, noise SD 2.09, top 30% = cases)."
- `simulation_log.csv`: per tag -> weight, causal SNP id/pos, own_causal, n_candidates, r2 with tag per ancestry.

## 7. SNP list, weights, annotations

GWAS p-values: [MISSING: no GWAS p-values anywhere — not in score.csv, raw PGS scoring file (columns: rsID, chr_name, effect_allele, other_allele, effect_weight, allelefrequency_effect, hm_source, hm_rsID, hm_chr, hm_pos, hm_inferOtherAllele), notebook, READMEs or report]

Nearest-gene annotations present in the repo (only these):
- rs7229639 -> "SMAD7" (REPORT §3.1). Weight 0.023, effect allele blank, other allele G, chr18:46450976 (GRCh37), published effect-allele freq 0.09.
- "FADS1/2 SNPs on chr11" (REPORT Q1; DREADME Limitations) — rsIDs not named there; the chr11 ~61.5 Mb tags are rs1535 (0.184), rs174537 (-0.095), rs174550 (0.134), rs4246215 (-0.165) per score.csv + simulation_log causal positions 11:61.57-61.59 Mb (mapping of these rsIDs to "FADS1/2" is by position inference, not stated in the repo).
- rs6983267 (MYC/CASC8/CCAT2, 8q24): weight 0.073, effect G / other T, chr8:128413305 GRCh37, published effect-allele freq 0.51; flipped in harmonisation (REF=G); simulation causal 8:128407443, r2 EUR 0.912 / AFR 0.377 / EAS 0.957 / SAS 0.882 / AMR 0.908 (n_candidates 2). [MISSING: no MYC/CASC8/CCAT2/8q24 gene annotation anywhere in the Lab 6 files]
- [MISSING: TCF7L2, BMP4 or any other nearest-gene annotation — not in repo]
- Report also mentions rs11064437 as largest weight (0.218), no gene named.

Full weights (`data/score.csv`; 62 rows; rs202110856 weight 0.02 only in raw file):

| variant | EA | OA | weight | | variant | EA | OA | weight |
|---|---|---|---|---|---|---|---|---|
| rs10911251 | A | C | 0.076 | | rs3802842 | C | A | 0.098 |
| rs6687758 | G | A | 0.056 | | rs10849432 | T | C | 0.069 |
| rs6691170 | T | G | 0.001 | | rs3217810 | T | C | 0.143 |
| rs72647484 | T | C | 0.001 | | rs10774214 | T | C | 0.040 |
| rs11903757 | C | T | 0.072 | | rs11064437 | C | T | 0.218 |
| rs812481 | G | C | 0.065 | | rs11169552 | C | T | 0.014 |
| rs35360328 | A | T | 0.127 | | rs7136702 | T | C | 0.034 |
| rs10936599 | C | T | 0.015 | | rs3184504 | C | T | 0.089 |
| rs647161 | A | C | 0.053 | | rs73208120 | G | T | 0.188 |
| rs1321311 | A | C | 0.078 | | rs1957636 | T | C | 0.047 |
| rs4711689 | A | G | 0.005 | | rs4444235 | C | T | 0.086 |
| rs7758229 | T | G | 0.039 | | rs4779584 | T | C | 0.105 |
| rs16892766 | C | A | 0.200 | | rs11632715 | A | G | 0.025 |
| rs2450115 | T | C | 0.042 | | rs16969681 | T | C | 0.046 |
| rs6469656 | A | G | -0.005 | | rs79900961 | G | A | 0.030 |
| rs10505477 | A | G | -0.008 | | rs9929218 | G | A | 0.080 |
| rs6983267 | G | T | 0.073 | | rs16941835 | C | G | 0.054 |
| rs7014346 | A | G | 0.096 | | rs12603526 | C | T | 0.068 |
| rs719725 | A | C | 0.035 | | rs4939827 | T | C | 0.111 |
| rs10904849 | G | T | 0.010 | | rs7229639 | (blank) | G | 0.023 |
| rs10795668 | G | A | 0.017 | | rs10411210 | C | T | 0.064 |
| rs704017 | G | A | 0.072 | | rs1800469 | G | A | 0.151 |
| rs1035209 | T | C | -0.011 | | rs2241714 | C | T | -0.124 |
| rs11190164 | G | A | 0.087 | | rs2423279 | C | T | 0.066 |
| rs4919687 | G | A | 0.043 | | rs4813802 | G | T | 0.076 |
| rs12241008 | C | T | 0.116 | | rs961253 | A | C | 0.094 |
| rs10506868 | T | C | -0.056 | | rs6066825 | A | G | 0.079 |
| rs11196172 | A | G | 0.062 | | rs4925386 | C | T | -0.001 |
| rs1535 | A | G | 0.184 | | rs6061231 | C | A | 0.042 |
| rs174537 | G | T | -0.095 | | | | | |
| rs174550 | T | C | 0.134 | | | | | |
| rs4246215 | G | T | -0.165 | | | | | |
| rs3824999 | G | T | 0.079 | | | | | |

## 8. Disagreements (listed, not reconciled)

1. `part3_crc/build.log` (gitignored) is from a different/earlier build and disagrees with the saved data, DREADME, NB and REPORT: "tags found in 1000 Genomes: 61 of 63 missing: ['rs7229639', 'rs202110856']" and "harmonization: keep 33, flip 28, drop 0" (vs 62 found, drop 1); "mean tag-causal r2: EUR 0.84, AFR 0.47, EAS 0.78, SAS 0.77, AMR 0.76" (vs 0.83/0.45/0.76/0.77/0.74); "cases: 0.300 overall; EUR 0.280, AFR 0.350, EAS 0.248, SAS 0.362, AMR 0.245" (vs meta.csv/NB 0.300/0.266/0.282/0.380/0.262).
2. REPORT Q3: "SAS people with low scores have about twice the risk the model gives them" (bottom two quintiles). NB cell 19: Q1 0.32 vs 0.15 (2.1x) but Q2 0.32 vs 0.23 (1.4x) — "about twice" holds for Q1 only.
3. NB markdown cells 0, 1 and 22 (professor's starter text, restored per commit 8fb87b2) still describe the Alzheimer's PGS000334 score, APOE caveat and "Explain why ... ε4"; the data/outputs are CRC PGS000148. Student note in cell 2 flags this.
4. Age: DREADME says "Age is also simulated (it is the Canvas cohort's age column)"; REPORT Setup says "The phenotype is simulated, and so is the age effect" (age effect vs age values); REPORT Q4 later says "simulated label and simulated age".
5. DREADME line "One row (none) publishes no `other_allele`" — template artifact: no row actually lacked other_allele (filled_other list empty), yet sentence asserts "One row".
6. REPORT table "% above pooled top-10% cutoff Overall 10.0" and "Mean tag-causal r²" row are not NB outputs (from definition and DREADME respectively; report states r² source).
7. TREADME/REPORT describe score as "63 SNPs"; NB proof block prints "score size : 62 variants" (explained in REPORT §3.1, consistent once explained).
8. REPORT calls the cutoff "top-10%"; the cutoff value 1.282 is an empirical 90th percentile, flagging 10.0% pooled — consistent.
No numeric disagreement found between REPORT tables and NB outputs for n, mean PRS, % flagged, % cases, AUC, O/E, intercept/slope, quintile values quoted, or CV AUCs.
