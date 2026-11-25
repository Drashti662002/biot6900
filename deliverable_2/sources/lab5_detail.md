# Lab 5 (Module 5, pharmacogenomics / VEPs): evidence inventory

Repo: `/Users/drashtibhanushali/biot6900/modules/module5/`. Nothing in the repo was modified. Recomputations used pandas on `data/*.tsv`, applying the notebook's own rules (`clinvar_label()`, `af_filter(MAX_AF=0.5)`, bands `[0,10,20,25,30,50]`, all in cell 0 / Q3 / Q4 of the graded notebook).

## 0. Provenance caveat (important)

- The **working-tree** `BIOT6900_Module5_Part3_Assignment.ipynb` has **no saved outputs**: all 9 code cells have `execution_count=None`. `git diff --stat` shows 573 lines deleted since HEAD, which is an uncommitted output strip.
- The **committed** version (`git show HEAD:./BIOT6900_Module5_Part3_Assignment.ipynb`, first added in commit ec47700) **does have outputs**. Every notebook number below comes from that HEAD copy. Cell indices are HEAD indices, which include markdown cells: Q1 = cell 5, Q2 = cell 7, Q3 = cell 9, Q4 = cell 11, Q5 = cell 13, Atlas = cell 15, discordance = cell 16.
- Part 2 copies: `BIOT6900_Module5_Part2_RealData (1).ipynb` has outputs in 15 of 39 cells. `...RealData.ipynb` and `...RealData (2).ipynb` are byte-identical (same md5) and have no outputs.
- `data/pgx_variants.tsv` (635 rows) and the per-gene TSVs hold the same DPYD variants, with 0 CADD mismatches. The `note` column marks the anchors.

## 1. CYP2C19 coverage (80 variants)

Source: HEAD notebook Q1 (cell 5) and Q2 (cell 7). Recomputed from `data/CYP2C19.tsv`, and the numbers match.

- Labels: other (drug response) 38, VUS 33, B/LB 9, P/LP 0.
- CADD 68/80 (85.0%), REVEL 54/80 (67.5%), AlphaMissense 54/80 (67.5%), gnomAD AF 55/80 (68.8%).

| class | n | CADD | REVEL | AM |
|---|---|---|---|---|
| missense | 54 | 53 | 53 | 54 |
| unannotated CNV / large deletion | 8 | 0 | 0 | 0 |
| synonymous | 7 | 6 | 0 | 0 |
| upstream | 5 | 4 | 0 | 0 |
| intergenic | 1 | 1 | 0 | 0 |
| intron | 1 | 1 | 0 | 0 |
| splice_donor | 1 | 0 | 0 | 0 |
| start_lost | 1 | 1 | **1** | 0 |
| stop_gained | 1 | 1 | 0 | 0 |
| stop_lost | 1 | 1 | 0 | 0 |

- Notebook output (Q2): "classes with ZERO REVEL scores (8)", "**25 variants** sit in classes no missense predictor can see." The report (section 2 table) says the same: 8 classes, 25 variants.
- Recomputed: classes with zero REVEL = 25 variants. Classes with zero AlphaMissense = **26**, because start_lost M1V (rs28399504) has REVEL 0.314 and no AM score. Rows with neither REVEL nor AM = 25. Rows with no tool at all = 11 of 80, and 10 after the AF filter (the `flag_no_tool_at_all` count in cell 16).

## 2. Pooled five-gene CADD calibration

Source: HEAD Q4 (cell 11). I recomputed it from `data/{CYP2C19,CYP2C9,VKORC1,DPYD,TPMT}.tsv` with MAX_AF 0.5 and got an **exact match** on every count.

Labelled variants with a CADD score, by gene:

| gene | B/LB | P/LP | total |
|---|---|---|---|
| CYP2C19 | 5 | 0 | 5 |
| CYP2C9 | 7 | 1 | 8 |
| DPYD | 91 | 39 | 130 |
| TPMT | 4 | 0 | 4 |
| VKORC1 | 11 | 9 | 20 |
| **All** | **118** | **49** | **167** |

| band | n | P/LP | B/LB | LR+ | notebook evidence | Haldane LR | report evidence |
|---|---|---|---|---|---|---|---|
| (0,10] | 76 | 2 | 74 | 0.07 | moderate, benign | 0.08 | moderate, benign |
| (10,20] | 34 | 1 | 33 | 0.07 | moderate, benign | 0.11 | moderate, benign |
| (20,25] | 17 | 8 | 9 | 2.14 | supporting, pathogenic | 2.14 | supporting, pathogenic |
| (25,30] | 6 | 4 | 2 | 4.82 | moderate, pathogenic (flag `<10 !`) | 4.31 | "under 10, not used" |
| (30,50] | 34 | 34 | 0 | inf | strong, pathogenic | 165.18 | strong, pathogenic |

No labelled variant falls outside the band range.

Per-gene contribution to each band (recomputed):

| band | CYP2C19 B | CYP2C9 B/P | DPYD B/P | TPMT B | VKORC1 B/P |
|---|---|---|---|---|---|
| (0,10] | 5 | 3/1 | 60/0 | 0 | 6/1 |
| (10,20] | 0 | 1/0 | 24/0 | 4 | 4/1 |
| (20,25] | 0 | 3/0 | 6/5 | 0 | 0/3 |
| (25,30] | 0 | 0/0 | 1/3 | 0 | 1/1 |
| (30,50] | 0 | 0/0 | 0/31 | 0 | 0/3 |

- The 3 low-band positives are not from DPYD. In (0,10]: CYP2C9 L208V rs72558191 (CADD 0.002) and VKORC1 −1639G>A rs9923231 (CADD 2.46, "Pathogenic; drug response"). In (10,20]: VKORC1 rs9934438 (CADD 10.42, "Pathogenic; drug response").
- Labelled variants dropped for lack of CADD: DPYD 66 P/LP + 6 B/LB; CYP2C19 2 B/LB; CYP2C9 1 B/LB; VKORC1 2 B/LB.

## 3. DPYD in detail

- `data/DPYD.tsv` has 426 rows. Six common intronic benign variants have AF above 0.5 and are removed, leaving 420. Labels: VUS 183, **P/LP 105**, B/LB 97, other 28, conflicting 6, none 1.
- Of the **105 P/LP**, **39 have a CADD score** (these are the "39 positives"), 7 have REVEL and 8 have AlphaMissense.
- **All 39 are recoverable at variant level**, each as a row with scores in both `data/DPYD.tsv` and `data/pgx_variants.tsv`. They are 39 distinct HGVS. 32 rows carry an rsID (31 distinct, because rs3918290 has two rows) and 7 have no rsID.
  - By consequence: 27 stop_gained, 7 missense, 3 splice_donor, 2 splice_acceptor.
  - By band: 31 in (30,50], 3 in (25,30], 5 in (20,25].
- The other 66 P/LP have no CADD: 27 frameshift, 13 splice_acceptor, 13 splice_donor, 11 with no consequence annotation, 1 stop_gained, 1 missense. They are still present as rows, but no calibration score exists for them.
- **Verdict:** DPYD is recoverable at variant level, not just as pooled band counts.

### Named DPYD alleles

| allele | rsID | present? | HGVS (hg38) | consequence | CADD | REVEL | AM | AF | ClinVar | stars | label → used in calibration? |
|---|---|---|---|---|---|---|---|---|---|---|---|
| *2A c.1905+1G>A | rs3918290 | yes, 2 rows | chr1:g.97450058C>T | splice_donor | 23.7 | – | – | 0.004555 | Pathogenic/Likely pathogenic | 2 | P/LP → yes, band (20,25] |
| (same rsID, other alt) | rs3918290 | | chr1:g.97450058C>G | splice_donor | 23.4 | – | – | 0.000008 | Likely pathogenic | 2 | P/LP → yes, band (20,25] |
| *13 c.1679T>G | rs55886062 | yes | chr1:g.97515787A>C (I560S) | missense | 31.0 | 0.920 | 0.8504 | 0.000388 | Pathogenic | 2 | P/LP → yes, band (30,50] |
| c.2846A>T | rs67376798 | yes | chr1:g.97082391T>A (D949V) | missense | 29.9 | 0.962 | 0.2925 | 0.00328 | Conflicting interpretations | 1 | conflicting → **excluded** |
| HapB3 c.1236G>A | rs56038477 | yes | chr1:g.97573863C>T | synonymous | 14.03 | – | – | 0.013112 | **Benign** | 2 | B/LB → **counted as a benign** in (10,20] |
| HapB3 c.1129-5923C>G | rs75017182 | yes | chr1:g.97579893G>C | intron | 0.378 | – | – | 0.013111 | drug response | 3 | other → excluded |
| (*9A c.85T>C, extra) | rs1801265 | yes | chr1:g.97883329A>G (C29R) | missense | NaN | NaN | 0.0426 | 0.273 | Pathogenic | 0 | P/LP but no CADD → excluded |

The C>T row is the *2A allele, since c.1905+1G>A on the minus strand corresponds to C>T on the plus strand.

### All 39 DPYD P/LP variants with CADD

Full table, also saved at `/private/tmp/claude-501/-Users-drashtibhanushali-biot6900-modules-module-10/ee468b5a-f628-4a63-9a5f-4d9be1868906/scratchpad/dpyd_pos39.csv`.

| hgvs_g | rsid | consequence | aa | CADD | band | REVEL | AM | AF | ClinVar | ★ |
|---|---|---|---|---|---|---|---|---|---|---|
| chr1:g.97515763C>T | rs2101975422 | stop_gained | W568X | 46.0 | 30-50 | | | | LP | 1 |
| chr1:g.97549582C>T | rs1214305633 | stop_gained | W501X | 44.0 | 30-50 | | | 4e-6 | LP | 1 |
| chr1:g.97515785G>A | rs1057516968 | stop_gained | R561X | 43.0 | 30-50 | | | 3.3e-5 | P/LP | 2 |
| chr1:g.97450101C>T | rs1057516388 | stop_gained | W621X | 42.0 | 30-50 | | | 7e-6 | LP | 2 |
| chr1:g.97549592G>A | – | stop_gained | Q498X | 42.0 | 30-50 | | | | LP | 1 |
| chr1:g.97573943C>A | rs78060119 | stop_gained | E386X | 42.0 | 30-50 | | | 7e-6 | P | 2 |
| chr1:g.97450133C>A | rs1057516671 | stop_gained | E611X | 41.0 | 30-50 | | | | LP | 1 |
| chr1:g.97549676G>A | – | stop_gained | Q470X | 41.0 | 30-50 | | | | LP | 1 |
| chr1:g.97699409G>A | rs1411354209 | stop_gained | R208X | 40.0 | 30-50 | | | 8e-6 | LP | 1 |
| chr1:g.97691782G>A | – | stop_gained | Q233X | 39.0 | 30-50 | | | | LP | 1 |
| chr1:g.97573856C>A | rs755031341 | stop_gained | E415X | 39.0 | 30-50 | | | 4e-6 | LP | 1 |
| chr1:g.97699370C>A | rs146170505 | stop_gained | E221X | 39.0 | 30-50 | | | 7e-6 | P/LP | 2 |
| chr1:g.97382432G>T | – | stop_gained | Y645X | 39.0 | 30-50 | | | | LP | 1 |
| chr1:g.97883353G>A | rs72549310 | stop_gained | R21X | 39.0 | 30-50 | | | 7e-6 | P/LP | 2 |
| chr1:g.97721529A>T | rs2101026231 | stop_gained | L155X | 38.0 | 30-50 | | | | P | 1 |
| chr1:g.97699523G>A | rs2100972370 | stop_gained | Q170X | 38.0 | 30-50 | | | | LP | 1 |
| chr1:g.97828127G>A | rs189768576 | stop_gained | R74X | 38.0 | 30-50 | | | 7.9e-5 | P/LP | 2 |
| chr1:g.97828142C>A | – | stop_gained | E69X | 38.0 | 30-50 | | | | P/LP | 2 |
| chr1:g.97699373G>A | rs1220608789 | stop_gained | Q220X | 37.0 | 30-50 | | | 4e-6 | LP | 1 |
| chr1:g.97828139G>A | rs141597515 | stop_gained | R70X | 37.0 | 30-50 | | | 2.6e-5 | P/LP | 2 |
| chr1:g.97828115T>A | rs776692894 | stop_gained | R78X | 37.0 | 30-50 | | | 4e-6 | LP | 1 |
| chr1:g.97082449G>A | rs2101558812 | stop_gained | Q930X | 36.0 | 30-50 | | | | LP | 1 |
| chr1:g.97595122G>A | rs2102259303 | stop_gained | Q299X | 35.0 | 30-50 | | | | LP | 1 |
| chr1:g.97593380G>T | rs768834079 | stop_gained | C322X | 35.0 | 30-50 | | | 7e-6 | LP | 1 |
| chr1:g.97679174G>T | rs754745863 | stop_gained | C257X | 35.0 | 30-50 | | | 1.3e-5 | LP | 1 |
| chr1:g.97515891G>C | – | stop_gained | Y525X | 35.0 | 30-50 | | | | LP | 1 |
| chr1:g.97883319G>T | – | stop_gained | S32X | 35.0 | 30-50 | | | | LP | 1 |
| chr1:g.97691776G>A | rs1801266 | missense | R235W | 34.0 | 30-50 | 0.952 | 0.8881 | 4.0e-5 | LP | 1 |
| chr1:g.97740456G>A | rs568132506 | missense | P86L | 34.0 | 30-50 | 0.846 | 0.8504 | 3.9e-5 | P/LP | 2 |
| chr1:g.97593274G>A | rs779557503 | missense | R358C | 33.0 | 30-50 | 0.747 | 0.2788 | 7e-6 | LP | 0 |
| chr1:g.97515787A>C | rs55886062 (*13) | missense | I560S | 31.0 | 30-50 | 0.920 | 0.8504 | 3.88e-4 | P | 2 |
| chr1:g.97699430T>G | rs72549308 | missense | S201R | 27.9 | 25-30 | 0.895 | 0.9982 | 7e-6 | LP | 1 |
| chr1:g.97740391C>T | rs746368304 | splice_donor | – | 27.9 | 25-30 | | | 1.3e-5 | P/LP | 2 |
| chr1:g.97082394A>G | rs372307932 | missense | I948T | 27.8 | 25-30 | 0.938 | 0.3934 | 4e-6 | LP | 1 |
| chr1:g.97079133T>A | rs72547602 | missense | D974V | 24.8 | 20-25 | 0.871 | 0.2727 | | P | 0 |
| chr1:g.97593389T>C | rs577302209 | splice_acceptor | – | 24.7 | 20-25 | | | 7e-6 | LP | 1 |
| chr1:g.97082471C>T | rs377221676 | splice_acceptor | – | 24.0 | 20-25 | | | 7e-6 | LP | 1 |
| chr1:g.97450058C>T | rs3918290 (*2A) | splice_donor | – | 23.7 | 20-25 | | | 4.555e-3 | P/LP | 2 |
| chr1:g.97450058C>G | rs3918290 | splice_donor | – | 23.4 | 20-25 | | | 8e-6 | LP | 2 |

The 66 P/LP rows without CADD are saved at `/private/tmp/claude-501/-Users-drashtibhanushali-biot6900-modules-module-10/ee468b5a-f628-4a63-9a5f-4d9be1868906/scratchpad/dpyd_plp_nocadd66.csv`.

The report itself has **no DPYD variant-level discussion**. It mentions DPYD only in the "39 of the 49 positives come from DPYD" bullet (section 1) and in the pooling sentence.

## 4. CYP2C19 anchors

Sources: `data/CYP2C19.tsv`, HEAD Q1 and Q3 outputs, and report section 3.

| allele | HGVS | consequence | CADD | band (LR) | REVEL | AM | AF | ClinVar | ★ |
|---|---|---|---|---|---|---|---|---|---|
| *2 rs4244285 | chr10:g.94781859G>A | synonymous | **5.686** | (0,10], LR 0.07 moderate benign | blank | blank | 0.167953 | drug response | 0 |
| *2 rs4244285, 2nd row | chr10:g.94781859G>T | synonymous | 2.606 | (0,10] | blank | blank | NaN | drug response | 0 |
| *17 rs12248560 | chr10:g.94761900C>T | upstream | **3.031** | (0,10] | blank | blank | 0.200120 | Benign; other | 2 |

- The report quotes CADD 5.69 / 5.7 for *2 and 3.03 / 3.0 for *17.
- *17 is one of the 5 CYP2C19 B/LB variants in the calibration set.
- *2 is labelled "other", so it is not in the calibration set. The report places it in the 0–10 band hypothetically.

## 5. VUS reclassification

Sources: HEAD Q5 (cell 13) and `pgx_classified_CYP2C19.csv` (33 rows), which agree with each other and with report and README.

- 33 VUS in total.
- 8 lean pathogenic, all in band (20,25] with LR 2.14 ("supporting").
- 16 lean benign: 8 in (0,10] and 8 in (10,20], LR 0.07.
- 9 unchanged: 3 in (25,30], which counts as insufficient data (n<10), and 6 with no CADD score.
- None is reclassified to P/LP or B/LB.

## 6. rs72558186 (*7)

- Data: chr10:g.94781999T>A, splice_donor_variant (`consequence_all` = intron_variant;splice_donor_variant), clinvar_id 39348, drug response, **4 stars**. CADD, REVEL, AM and gnomAD are all NaN.
- Report section 2: "my unscored variant", exon 5/intron 5 donor, trips `flag_no_tool_at_all`. Proposed next steps are SpliceAI donor-loss (>0.8) and then a minigene assay or carrier-liver RNA-seq. The skipped Atlas could have given a MERGED_SPLICING attribution.
- The report's claim that every column is empty matches the data.

## 7. UGT1A1 / irinotecan / fluoropyrimidine / CPIC

- **UGT1A1, irinotecan, 5-FU, capecitabine, DPWG: no mention** anywhere in the student's files: report (md and docx), README, Part 3 notebook, Part 2 notebooks, or Part 1 notebook.
- **Fluoropyrimidine:** appears only in the course-provided `BIOT6900_Module5_Part2_LabGuide.pdf`, in its gene table ("DPYD fluoropyrimidines … rare but severe toxicity"). The student never writes it.
- **CPIC:** the report (section 3) mentions it once, about *2: "CPIC practice-guideline status and drives a clopidogrel boxed warning". The Part 2 notebook's warfarin section mentions CPIC's warfarin algorithm, which is course text.

## 8. Disagreements and inconsistencies (listed, not reconciled)

1. **Outputs stripped from the graded notebook.** The working tree has 0 outputs. The outputs that support every reported number exist only in git HEAD.
2. **Two different "25 non-missense" sets.** Q2 (cell 7) counts 25 among the 80 variants, using classes with zero REVEL: it excludes start_lost and includes 7 synonymous. The Atlas cell (15) counts 25 among the 78 filtered variants: it includes start_lost, has 6 synonymous, and includes 8 CNVs not shown in its value_counts. The report treats them as one set ("the 25 non-missense variants above are exactly what the Atlas is for").
3. **"No missense predictor scored" holds for REVEL only.** AlphaMissense is blank for 26 variants (start_lost M1V has REVEL 0.314 but no AM). The report frames the 25 as classes "no missense predictor scored" / "REVEL returned nothing".
4. **(25,30] band label.** The notebook labels it "moderate, pathogenic" with a `<10 !` flag. The report table says "under 10, not used". The underlying numbers are the same.
5. ***2 second row ignored.** *2 has a second row (G>T, CADD 2.606, no AF). The report gives only 5.69, though it does mention the two rows.
6. **"39 of 49 positives DPYD" is CADD-conditional.** DPYD has 105 P/LP after filtering, and 66 of them (mostly frameshift and splice) are silently dropped for lack of CADD. The report does not mention this.
7. **HapB3 c.1236G>A counted benign; c.2846A>T excluded.** HapB3 c.1236G>A (rs56038477) is ClinVar "Benign" and enters the pooled calibration as a benign in band (10,20]. DPYD c.2846A>T, a lab anchor, is "Conflicting" and excluded. The report does not discuss either.
8. **Part 2 calibration differs.** The Part 2 (1) notebook's CADD calibration (P/LP 49, B/LB 125; LRs 0.06 / 0.08 / 2.27 / 5.10 / inf) differs from the report's (118 B/LB). It comes from a different notebook and filter, so it is not a direct contradiction, but these are different numbers for "the same" calibration.
9. **Calibration counts.** No disagreement: the report's table, 167/49/118, VUS counts, filter counts (2 removed; 19/14 for the 0.001 filter) and coverage percentages all match the HEAD notebook outputs and my recomputation from the data.

## MISSING

- [MISSING: DPYD HapB3 or c.2846A>T discussion in the report]
- [MISSING: any UGT1A1 / irinotecan content]
- [MISSING: executed outputs in the working-tree Part 3 notebook; they are available only in git HEAD]
