# Portfolio Deliverable 2 — Phase 1 Inventory

Colon adenocarcinoma (COAD) case study. Source material: Labs 5, 6, 7, 8 and Portfolio Deliverable 1.
Compiled 2026-11-24. Nothing in the lab folders was modified or re-executed.

Per-lab detail, with every value and its source location, is in `sources/`:
`lab5_detail.md`, `lab6_detail.md`, `lab7_deliverable1_detail.md`, `lab8_detail.md`,
plus `dpyd_pos39.csv` and `dpyd_plp_nocadd66.csv`.

---

## 0. Stop conditions

| # | Condition | Status |
|---|---|---|
| 1 | DPYD variant-level data recoverable | **CLEARED.** Full per-variant rows with scores in `modules/module5/data/DPYD.tsv`. Verified directly. |
| 2 | Per-ancestry calibration quintiles found | **CLEARED.** Printed as text, Lab 6 notebook cell 19. Figure 2B can be redrawn exactly. |
| 3 | Lab 7 graph output and D1 target list found | **CLEARED.** 54 nodes / 61 edges in CSV; full PageRank and community output; D1 15-gene triage table present. |

Three findings need a decision before Section 2 and Section 5 are written. See §7.

**Repo protection (done before anything else).** The Lab 5 graded notebook had an uncommitted
working-tree change that had deleted all saved output: 0 of 9 code cells held outputs on disk
(`git diff --stat`: 34 insertions, 573 deletions). That change is now stashed as
`stash@{0}` = `8b74674a0915607922c5883868c931eae4008323`
("wip: lab5 notebook outputs cleared"). The notebook on disk is now byte-identical to HEAD and all
9 code cells hold outputs, execution counts 1–9 in order. The notebook was never executed or
re-saved; it was read as a file only. Recover the stashed change with
`git stash pop` or inspect it with `git stash show -p stash@{0}`.

---

## 1. Lab 5 — Variant interpretation and pharmacogenomics

Source files: `modules/module5/data/{DPYD,CYP2C19,CYP2C9,TPMT,VKORC1,pgx_variants}.tsv`,
`pgx_classified_CYP2C19.csv`, `Module5_Assignment_Report.md`,
`BIOT6900_Module5_Part3_Assignment.ipynb`.

**Notebook output location.** The working-tree copy of the graded notebook has no saved outputs:
an uncommitted change deleted 573 lines of output. Every notebook number below is read from the
committed copy, `git show HEAD:"modules/module5/BIOT6900_Module5_Part3_Assignment.ipynb"`.

### 1.1 CYP2C19 tool coverage (80 ClinVar variants)

| Tool | Variants scored |
|---|---|
| CADD | 68 |
| REVEL | 54 |
| AlphaMissense | 54 |
| gnomAD frequency | 55 |

25 variants across 8 consequence classes have no REVEL score (committed cell 7; report §2).
AlphaMissense is blank for 26, not 25: the start-lost variant M1V has REVEL 0.314 but no
AlphaMissense score. The report's phrase "no missense predictor scored" is exact for REVEL only.

### 1.2 Pooled five-gene CADD calibration

167 labelled variants, 49 P/LP and 118 B/LB. Recomputed from the five gene files; matches the
report exactly.

| CADD band | n | P/LP | B/LB | LR+ |
|---|---|---|---|---|
| 0–10 | 76 | 2 | 74 | 0.07 |
| 10–20 | 34 | 1 | 33 | 0.07 |
| 20–25 | 17 | 8 | 9 | 2.14 |
| 25–30 | 6 | 4 | 2 | 4.82 (n < 10) |
| 30–50 | 34 | 34 | 0 | inf (165.18 with Haldane correction) |

Per gene: DPYD 130 labelled (39 P/LP, 91 B/LB); VKORC1 20 (9/11); CYP2C9 8 (1/7);
CYP2C19 5 (0/5); TPMT 4 (0/4). DPYD supplies 31 of the 34 positives in the 30–50 band.
None of the three positives below CADD 20 is a DPYD variant: they are CYP2C9 L208V,
VKORC1 rs9923231 and VKORC1 rs9934438.

### 1.3 DPYD at variant level — RECOVERABLE

`data/DPYD.tsv` holds 426 rows with `rsid`, `hgvs_g`, `consequence`, `cadd_phred`, `revel`,
`alphamissense`, `gnomad_af`, `clinvar_sig`, `clinvar_stars`. All 39 P/LP variants counted in the
calibration are present as full rows: 27 stop-gained, 7 missense, 3 splice-donor, 2 splice-acceptor;
7 carry no rsID. Listed in `sources/dpyd_pos39.csv`.

Star alleles, verified directly against the file:

| Allele | rsID | Consequence | CADD | REVEL | AlphaMissense | ClinVar | In calibration? |
|---|---|---|---|---|---|---|---|
| *2A | rs3918290 (C>T) | splice donor | 23.7 | — | — | P/LP, 2★ | yes, band 20–25 |
| *2A second row | rs3918290 (C>G) | splice donor | 23.4 | — | — | Likely pathogenic, 2★ | yes |
| *13 | rs55886062 | missense I560S | 31.0 | 0.92 | 0.850 | Pathogenic, 2★ | yes, band 30–50 |
| c.2846A>T | rs67376798 | missense | 29.9 | 0.962 | 0.293 | Conflicting, 1★ | **no — excluded** |
| HapB3 c.1236G>A | rs56038477 | synonymous | 14.03 | — | — | Benign, 2★ | **yes, as a benign** |
| HapB3 intronic | rs75017182 | intron | 0.378 | — | — | drug response, 3★ | no |

### 1.3a The 39 is a filtered subset — named finding

`data/DPYD.tsv` holds **105** P/LP DPYD variants. Only **39 carry a CADD score**, so the
calibration silently drops **66 of 105 (63%)** (`sources/dpyd_plp_nocadd66.csv`).

| Consequence | Dropped (no CADD), n = 66 | Scored (in calibration), n = 39 |
|---|---|---|
| frameshift | 27 | 0 |
| splice acceptor | 13 | 2 |
| splice donor | 13 | 3 |
| stop gained | 1 | 27 |
| missense | 1 | 7 |

By consequence class, CADD scored **27 of 28 stop-gained** DPYD positives, **2 of 40 splice**
variants and **0 of 27 frameshifts**. The failure is representation, not ranking: frameshifts are
indels, and a scoring scheme built over single-nucleotide variants has nothing to assign them.

Composition of the two sets:

- **39 scored positives:** 69% stop-gained (27), 18% missense (7), 13% splice (5)
- **66 dropped positives:** 41% frameshift (27), 39% splice (26), remainder missense and stop-gained

**Therefore the 30–50 band does not show that CADD separates pathogenic from benign in DPYD.**
That band holds 34 P/LP against 0 benign and returns an infinite LR+. What it shows is that CADD can
recognise a premature stop codon, which is the case requiring no model at all: the consequence
annotation already carries that information. The variants a clinic actually struggles with — 40
splice and 27 frameshift — are absent from the table that licenses the likelihood ratio.

This makes the infinite LR a **stronger problem than the Module 5 report stated**. The report treats
the empty benign cell as a small-sample artifact to be corrected with a Haldane adjustment
(165.18). The real issue is that the empty cell sits in a band populated by the easiest positives,
with the hard ones excluded by file format rather than by any analytic choice.

**c.2846A>T as label noise.** rs67376798 is an established 5-FU toxicity variant (CADD 29.9,
REVEL 0.962, AlphaMissense 0.293) excluded from the positive set because ClinVar lists it
"Conflicting interpretations of pathogenicity" (1★). Together with the 66, this is the second way
the denominator was chosen by the data rather than by the analyst: once by file format, once by
label. Note also the predictor disagreement on this variant, with CADD and REVEL high while
AlphaMissense reads 0.293.

### 1.3b HapB3 — the tag and the causal variant, both rows verbatim

| Field | c.1236G>A (exonic tag) | c.1129-5923C>G (intronic, causal) |
|---|---|---|
| rsID | rs56038477 | rs75017182 |
| hgvs_g | chr1:g.97573863C>T | chr1:g.97579893G>C |
| consequence | synonymous_variant | intron_variant |
| CADD | 14.03 | 0.378 |
| REVEL | — | — |
| AlphaMissense | — | — |
| ClinVar | **Benign**, 2★ | drug response, 3★ |
| In calibration | **yes, as a benign variant in band 10–20** | no |

The exonic variant is a tag in linkage disequilibrium with the intronic variant that does the
damage. The intronic variant creates a cryptic splice site and causes exon 11 skipping; CPIC
assigns HapB3 a reduced-function activity-score penalty. The calibration scores the tag, calls it
benign on ClinVar's label, and never sees the causal allele, which sits at CADD 0.378 because
intronic positions under weak constraint are exactly what CADD ranks low.

### 1.3c c.2846A>T — excluded despite being established

rs67376798, CADD 29.9, REVEL 0.962, AlphaMissense 0.293, ClinVar "Conflicting interpretations of
pathogenicity" (1★). It is excluded from the positive set on the ClinVar label alone, although it
is an established 5-FU toxicity variant with a CPIC activity-score penalty. Note the internal
disagreement among the predictors: CADD and REVEL both score it high while AlphaMissense scores it
0.293, which is low.

The Lab 5 report says nothing about DPYD beyond "39 of the 49 positives come from DPYD". There is
no star-allele, activity-score or HapB3 discussion anywhere in the student's work.

### 1.4 CYP2C19 anchors

| Allele | rsID | CADD | Band | LR+ | REVEL | AlphaMissense | ClinVar |
|---|---|---|---|---|---|---|---|
| *2 | rs4244285 | 5.686 | 0–10 | 0.07 | — | — | drug response, 0★ |
| *17 | rs12248560 | 3.031 | 0–10 | 0.07 | — | — | Benign |

*2 is synonymous and creates a cryptic splice donor. A second *2 row (G>T, CADD 2.606) exists and
is not used in the report. *17 is one of CYP2C19's five benign calibration variants.

### 1.5 VUS reclassification

33 VUS total: 8 lean pathogenic (band 20–25), 16 lean benign, 9 unchanged (3 in band 25–30,
6 unscored). Notebook, `pgx_classified_CYP2C19.csv`, report and README all agree.

### 1.6 rs72558186 (*7)

Splice-donor variant, 4★, drug response. CADD, REVEL, AlphaMissense and gnomAD all empty, so it
trips the `flag_no_tool_at_all` check. The report proposes SpliceAI, then a minigene assay or RNA-seq.

### 1.7 Coverage gaps for Section 2

UGT1A1, irinotecan, 5-FU and capecitabine appear nowhere in the student's Lab 5 work.
"Fluoropyrimidines" appears only in the course lab-guide PDF. CPIC appears once in the report,
about *2 and clopidogrel. The CPIC DPYD recommendation and the UGT1A1 gap both need external
citations; neither can be sourced from this lab.

---

## 2. Lab 6 — Polygenic risk

Source files: `modules/module 6/part3_crc/data/{score,variants,meta,genotypes,simulation_log}.csv`,
`data/README.md`, `BIOT6900_Module6_Part3_Graded_CRC.ipynb`, `module 6_report.docx`.

### 2.1 What is real and what is simulated

Real: genotypes from 1000 Genomes phase 3 (GRCh37); weights and alleles from PGS Catalog
PGS000148 (Jeon 2018, log-OR from a European-only GWAS); 2,336 people with ancestry labels.
**Simulated: case status and age.** The data README states "THE PHENOTYPE IS SIMULATED … Age is
also simulated." Recipe, seed 6900: each tag gets a causal SNP within ±250 kb at r² > 0.6 in EUR;
liability = genetic score (scaled to SD 1) + 0.35 × (age − 60) / 8 + N(0, 2.09²); the top 30% of
the cohort are cases.

Every AUC, calibration figure and O/E below therefore depends on a simulated label.

### 2.2 Harmonisation (cell 10)

61 variants used, 28 allele-flipped, 1 dropped. `score.csv` holds 62 rows.
rs7229639 (SMAD7 locus) was dropped because its published effect allele is blank — verified:
the row exists with `effect_allele` NaN, `other_allele` G, weight 0.023. It is still included in
the simulated disease model. rs202110856 has no 1000 Genomes phase 3 record, so it was never in
the data.

### 2.3 Per-ancestry performance

Cutoff PRS > 1.282 is the empirical 90th percentile of the z-scored PRS across the whole cohort
(cell 14). n, mean PRS, flag rate and case rate from cell 14; AUC from cell 16; O/E from cell 19;
mean tag-causal r² from `data/README.md` (not printed in the notebook; column means of
`simulation_log.csv` agree).

| Ancestry | n | Mean PRS | % flagged | % cases | AUC | Mean tag-causal r² | O/E |
|---|---|---|---|---|---|---|---|
| EUR | 500 | +0.16 | 12.2 | 30.0 | 0.685 | 0.83 | 1.00 |
| AFR | 500 | −0.09 | 8.4 | 26.6 | 0.567 | 0.45 | 0.99 |
| EAS | 500 | −0.17 | 7.2 | 28.2 | 0.689 | 0.76 | 1.10 |
| SAS | 489 | +0.30 | 14.9 | 38.0 | 0.606 | 0.77 | 1.20 |
| AMR | 347 | −0.28 | 6.3 | 26.2 | 0.670 | 0.74 | 1.08 |

Overall AUC 0.647. EUR-fitted logistic model: intercept −1.054, slope 0.720 (cell 19).

### 2.4 Calibration quintiles (cell 19) — Figure 2B source

Quintiles taken within each ancestry on risk predicted by the EUR-fitted model.

| | Q1 | Q2 | Q3 | Q4 | Q5 |
|---|---|---|---|---|---|
| EUR predicted | .13 | .22 | .29 | .36 | .51 |
| EUR observed | .13 | .23 | .26 | .36 | .52 |
| AFR predicted | .11 | .18 | .25 | .32 | .47 |
| AFR observed | .16 | .26 | .31 | .26 | .34 |
| EAS predicted | .11 | .17 | .23 | .31 | .45 |
| EAS observed | .10 | .24 | .23 | .34 | .50 |
| SAS predicted | .15 | .23 | .30 | .38 | .53 |
| SAS observed | .32 | .32 | .30 | .39 | .58 |
| AMR predicted | .10 | .16 | .22 | .29 | .45 |
| AMR observed | .11 | .20 | .25 | .29 | .46 |

AFR observed rises .16 → .34 while predicted rises .11 → .47: the curve is flat, the model barely
ranks anyone, yet average calibration is O/E 0.99.

**SAS under-prediction is one quintile, not two.** Q1 observed .32 vs predicted .15 is 2.1×.
Q2 is .32 vs .23, which is 1.4×. The brief's "bottom quintiles … roughly two-fold" holds for Q1 only.

### 2.5 Incremental value over age (cell 21)

5-fold cross-validated AUC: age alone 0.579, age + PRS 0.666. Both the age effect and the outcome
are simulated.

### 2.6 SNP annotation gap

`score.csv` carries `variant`, `effect_allele`, `other_allele`, `weight` and nothing else.
**No GWAS p-values exist anywhere in Lab 6** — not in the score file, the notebook, the READMEs or
the report. No gene labels either, except SMAD7 for rs7229639 and "FADS1/2 SNPs on chr11" with no
rsIDs. rs6983267 is present with weight 0.073 (effect allele G, other T, flipped) but carries no
MYC/CASC8/CCAT2 annotation. All gene-to-SNP mapping for Section 4 must come from Lab 7's
`gwas_lookup.csv`.

Weights of the four SNPs that Lab 7 also carries as variant nodes, verified in `score.csv`:
rs6983267 0.073, rs4939827 0.111, rs4444235 0.086, rs11196172 0.062.

---

## 3. Lab 7 — Knowledge graph

Source files: `modules/module 7/{nodes.csv,edges.csv,data/gwas_lookup.csv,data/targets_coad.csv}`,
`Module7_report.pdf`, `BIOT6900_Module7_Starter.ipynb`.

**Notebook output location.** The working-tree notebook is uncommitted and partly un-run: cells
3 and 7 fail without a database, cells 10–18 have no outputs, and it lacks the `CONFIDENCE` dict,
so it cannot have produced the submitted `edges.csv`. All executed output below comes from the
committed copy at `bf72daa`, which is fully run and matches the CSVs and the report.

### 3.1 Nodes and edges

54 nodes, 61 edges, verified in the CSVs and matching the report.
Nodes: 18 gene (15 COAD plus IL23R, JAK2, IL12B), 2 disease, 5 phenotype, 9 pathway, 13 drug,
7 variant.

Edges by source and confidence (`confidence` comes from a dict keyed on `source`):

| Source | n | Confidence |
|---|---|---|
| GWAS (COAD) | 4 | high |
| GWAS Catalog (UC) | 3 | high |
| CPTAC somatic | 2 | high |
| multi-omics (TCGA/CPTAC) | 4 | medium |
| pharmacology | 5 | **low** |
| HPO | 8 | low |
| Reactome | 15 | low |
| ChEMBL | 13 | high |
| GWAS Catalog (Module 6 PRS) OCCURS_IN | 4 | low |
| CPTAC somatic OCCURS_IN | 3 | high |

Totals: 25 high, 4 medium, 32 low.

### 3.2 PageRank

KRAS 0.0377; BRAF, EGFR, TYMS 0.0231; ten genes at 0.0158; DES, MYH11, IL23R, CR2 at the 0.0086
floor. Genes gain PageRank only from inbound drug and variant edges, so the ranking partly measures
how many drugs and variants were attached. MYC carries the strongest genetic evidence in the graph
(p = 2 × 10⁻¹²²) at mid-table PageRank. The report's Table 2 drops JAK2 and IL23R from the output.

### 3.3 Communities — all 10, full membership

`greedy_modularity_communities`, committed notebook cell 16. Part I targets in **bold**.

| # | Members |
|---|---|
| 0 | Abdominal pain, Hematochezia, IL12B, IL23R, Iron deficiency anemia, JAK2, Tofacitinib, Ulcerative colitis, Ustekinumab |
| **1** | **CR2**, Colon adenocarcinoma, Complement cascade, Constipation, **DES**, **MYH11**, Smooth muscle contraction, Weight loss |
| 2 | BRAF, BRAF V600E, Cetuximab, EGFR, EGFR-RAS-MAPK signaling, Encorafenib, Panitumumab |
| 3 | Capecitabine, Fluorouracil, Irinotecan, Nucleotide synthesis & DNA replication, TOP1, TYMS |
| 4 | Adagrasib, KRAS, KRAS G12C, KRAS G12D, Sotorasib |
| 5 | MYC, TCF7L2, WNT/beta-catenin signaling, rs11196172, rs6983267 |
| 6 | BMP4, SMAD7, TGF-beta/BMP signaling, rs4444235, rs4939827 |
| 7 | Bevacizumab, VEGF signaling / angiogenesis, VEGFA |
| 8 | Aspirin, PTGS2, Prostaglandin synthesis |
| **9** | Acetazolamide, **CA2**, Carbon dioxide hydration |

**Community 1 is the one the report reads as tissue composition**, calling it "a warning about my own
Module 2 ranking". Two separate points, to be kept distinct:

1. **The gene membership** is CR2, DES and MYH11, three down-in-tumour genes whose proteins are
   smooth muscle (DES, MYH11) and a B-cell/FDC receptor (CR2). Tissue composition is the most likely
   interpretation of that membership.
2. **The cluster also contains** the COAD disease node and two phenotypes, Constipation and Weight
   loss, which arrive on hand-entered HPO edges marked low confidence. That is a weakness in how the
   cluster was assembled and is not evidence about tissue composition.

**CA2 is in community 9, alone with acetazolamide and carbon dioxide hydration.** Verified in
`edges.csv`: CA2's only pathway edge is Carbon dioxide hydration. The community result supports a
tissue-composition reading for **CR2, DES and MYH11 only**. CA2 must never be written into that
argument; its weakness is direction, not cell type. See §7.4.

Three of the four Part I genes carried into Part II therefore sit in one community, and the fourth,
the Part I top-ranked gene, sits in its own.

### 3.3a The 2-hop drug result — all 11 rows

Committed notebook cell 10. Query:
`MATCH (d:Drug)-[:TARGETS]->(g:Gene)-[:ASSOCIATED_WITH]->(dis:Disease) WHERE dis.name = $disease`

| Drug | Gene | Status |
|---|---|---|
| Acetazolamide | **CA2** | approved (other indication) |
| Adagrasib | KRAS | approved (CRC) |
| Aspirin | PTGS2 | approved (other indication) |
| Bevacizumab | VEGFA | approved (CRC) |
| Capecitabine | TYMS | approved (CRC) |
| Cetuximab | EGFR | approved (CRC) |
| Encorafenib | BRAF | approved (CRC) |
| Fluorouracil | TYMS | approved (CRC) |
| Irinotecan | TOP1 | approved (CRC) |
| Panitumumab | EGFR | approved (CRC) |
| Sotorasib | KRAS | approved (CRC) |

Nine of eleven are approved in CRC. The two "other indication" rows are the repurposing-shaped ones,
and one of them is CA2 with acetazolamide. **CR2, DES and MYH11 return no drug at all** — they have
no TARGETS edge, which is why they sit at the PageRank floor.

The separate 5-hop repurposing query (cell 11) returned exactly two candidates, both via ulcerative
colitis and the three shared symptoms: Tofacitinib (JAK2) and Ustekinumab (IL12B).

### 3.4 Drug–target edges behind the clinical claims

| Drug | Target | Source | Confidence |
|---|---|---|---|
| Fluorouracil | TYMS | ChEMBL | high |
| Capecitabine | TYMS | ChEMBL | high |
| Irinotecan | TOP1 | ChEMBL | high |
| Acetazolamide | CA2 | ChEMBL | high |

The gene–disease edges for TYMS and TOP1 are `pharmacology`, confidence **low**: the gene is in the
graph because a COAD drug targets it, so a path from drug to gene to disease is circular. The report
says so. The 2-hop drug query returned 11 drug–gene rows.

### 3.5 The IL12B path

Captured edge by edge with source and confidence: ChEMBL high → GWAS Catalog high (p 1 × 10⁻²¹) →
HPO low → HPO low. The two HPO phenotype edges carrying the repurposing claim were hand-entered and
are low confidence. The p40-versus-p19 complication is discussed, mirikizumab is named as the
cleaner test, and an AOM/DSS three-arm plan is proposed.

### 3.6 Cross-lab variant nodes

| Variant | Edge | Confidence | Also in Lab 6? | Lab 6 weight |
|---|---|---|---|---|
| rs6983267 | OCCURS_IN → MYC | low | yes | 0.073 |
| rs4939827 | OCCURS_IN → SMAD7 | low | yes | 0.111 |
| rs4444235 | OCCURS_IN → BMP4 | low | yes | 0.086 |
| rs11196172 | OCCURS_IN → TCF7L2 | low | yes | 0.062 |

### 3.6a The MYC p-value — FOUND, fully sourced

Print convention: the document prints SMAD7's p-value as 3 × 10⁻³⁰. The source file's own float
artifact (`3.0000000000000003e-30`) stays verbatim in the data files, `evidence_table.csv` and the
provenance log.

`modules/module 7/data/gwas_lookup.csv`, row 1, all columns verbatim:

| Column | Value |
|---|---|
| gene | MYC |
| gwas_snp | rs6983267 |
| gwas_trait | colorectal cancer |
| gwas_p | **2e-122** |
| gwas_pmid | 36539618 |
| gwas_study | GCST90129505 |
| gwas_author | Fernandez-Rozadilla C |
| gwas_mapped_genes | **CASC8;CCAT2;PCAT1;POU5F1B** |

The same value appears on the MYC → Colon adenocarcinoma edge in `edges.csv`
(`pvalue` 2.000000e-122, `confidence` high, `source` GWAS). Not `[MISSING]`.

### 3.6b rs6983267 — both mappings recorded precisely

| Lab | What the file holds | Source |
|---|---|---|
| Lab 6 | `variant` rs6983267, `effect_allele` G, `other_allele` T, `weight` 0.073, allele-flipped during harmonisation. **No gene column exists in the file.** No MYC, CASC8, CCAT2 or 8q24 label anywhere in Lab 6. | `module 6/part3_crc/data/score.csv` |
| Lab 7, `gene` column | **MYC** | `module 7/data/gwas_lookup.csv` |
| Lab 7, `gwas_mapped_genes` column | **CASC8;CCAT2;PCAT1;POU5F1B** — MYC is not among them | same file, same row |
| Lab 7, graph edge | `rs6983267 OCCURS_IN MYC`, confidence **low**, source "GWAS Catalog (Module 6 PRS)" | `module 7/edges.csv` |

The disagreement is inside one file and one row: the `gene` column assigns the SNP to MYC while the
Catalog's `gwas_mapped_genes` column, in the adjacent cell, assigns it to four non-coding RNA genes
and a pseudogene, with MYC absent from that list. The same pattern appears on rs4444235: `gene` says
BMP4, `gwas_mapped_genes` says DDHD1-DT;MIR5580. For rs4939827 (SMAD7) and rs11196172 (TCF7L2) the
two columns agree.

**The graph recorded its own uncertainty correctly.** The `rs6983267 OCCURS_IN MYC` edge is marked
confidence **low**, which is the right label for an assignment the Catalog does not support. The
error is not that the graph hid the judgement call; it is that the MYC node, which carries the
strongest genetic association in the whole graph at p = 2 × 10⁻¹²², exists on that call.

**The cross-lab asymmetry is the finding.** Lab 6's `score.csv` has no gene column at all: columns
are `variant`, `effect_allele`, `other_allele`, `weight` and nothing else. A PRS therefore scores a
patient on rs6983267 without ever committing to which gene the variant acts through, because it does
not need to. Lab 7 had to commit, because a graph node must be a named entity, and it committed to a
gene the Catalog does not list. Two labs, the same unknown, two different ways of handling it: one
declines the question, the other answers it by assumption and flags the answer as low confidence.

**rs7229639 is absent from Lab 7 entirely.** It exists only in Lab 6. The brief's planned sentence
needs care: SMAD7 is still a graph node and still carries its own PRS tag SNP (rs4939827, weight
0.111) plus a high-confidence GWAS edge at p = 3 × 10⁻³⁰. The dropped variant did not weaken the
SMAD7 node; it removed one of two tags for that locus. See §7.3.

### 3.7 Figure source

`my_knowledge_graph.png` came from committed cell 18, a layered LAYER-based layout with a legend.
The working copy's spring layout is different code. `nodes.csv` and `edges.csv` are enough to
redraw topology, colour and size; the layout code exists only in git.

---

## 4. Lab 8 — Cell-type resolution

Source files: `modules/module 8/BIOT6900_Module8_Assignment_Starter.ipynb` (NB, 0-based cells),
`Module8_ERBB2.pdf`, `Module8_report.docx`.

### 4.1 Datasets

Single cell: CZ CELLxGENE Census, stable release 2025-11-08, filtered to breast cancer, primary
data only. 1,101,430 cells matched; 10 donors × 2,000 cells sampled; 18,388 cells after QC
(NB 5, 7, 10). Assays: 10x 3' v3 12,635; 10x 3' v2 3,753; 10x 5' v2 2,000.
Donor IDs appear only in an image legend (NB 10): HTAPP-211, -313, -516, -806, -889, -997,
QianBRCA_BC_5, pal_Patient 0056, 0176, 0178. There is no named atlas or collection — it is a
pooled Census query.

Spatial: 10x Visium V1_Breast_Cancer_Block_A_Section_1, 3,798 spots × 36,601 genes; the 500-count
spot filter removed none (NB 26, 30).

**The disease is breast cancer, not COAD.** This is the method-demonstration section.

### 4.2 Integration (NB 18)

Silhouette scores on a 5,000-cell subsample. No other metric was computed.

| Method | Batch silhouette (lower better) | Cell-type silhouette (higher better) |
|---|---|---|
| None | 0.095 | −0.016 |
| Harmony | −0.082 | 0.006 |
| scVI | −0.028 | 0.035 |

scVI chosen, 60 epochs on CPU.

### 4.3 Label correction (NB 10)

"malignant cell" contains the substring "t cell", so the provided labelling function assigned all
8,526 malignant cells to T. They were moved to Epithelial/Tumor. After correction: Tumor 8,678,
T 2,354, Other 2,976.

### 4.4 ERBB2 by cell type (NB 23) — OUTPUT-BACKED, full table

Printed in full by cell 23 as a text output, so Figure 4A can be redrawn exactly from numbers.
Computed on `adata.raw` (all genes, log-normalised) grouped by the corrected `cell_type_broad`.
The UMAP in the same cell is image-only.

| Cell type | n cells | Mean log-expression | Fraction expressing |
|---|---|---|---|
| **Epithelial/Tumor** | **8,678** | **0.633** | **0.364** |
| NK | 50 | 0.258 | 0.100 |
| DC | 36 | 0.217 | 0.139 |
| Other | 2,976 | 0.192 | 0.171 |
| B | 529 | 0.134 | 0.072 |
| Fibroblast | 1,072 | 0.112 | 0.096 |
| T | 2,354 | 0.056 | 0.025 |
| Macrophage | 1,189 | 0.048 | 0.030 |
| Endothelial | 1,431 | 0.042 | 0.038 |
| Monocyte | 73 | 0.000 | 0.000 |

Section 5's core contrast is therefore fully sourced and does not need rebuilding.

**One precision note for the write-up.** The Lab 8 report says "every other type is at or below 0.26,
and only 2.5% of T cells, 3% of macrophages and 4% of endothelial cells show it". The 0.26 ceiling
is correct (NK 0.258 is the highest non-tumour mean). But the 2.5–4% range covers only the three
types named: DC is 13.9% and "Other" is 17.1%. Both are small or heterogeneous groups (36 and 2,976
cells), so the argument holds, but the write-up should say "the named immune and endothelial types"
rather than implying every non-tumour type sits at 2.5–4%.

### 4.5 Spatial

Leiden resolution 0.6 gave 9 domains (NB 28). Domains 2–8 are tumour nests with ERBB2 means
1.31–1.65; domains 0–1 are stroma and immune (NB 30). Per domain (NB 39), 96–100% of spots in
nests 2–8 express ERBB2, against 80% in domain 0 (mean 1.22) and 84% in domain 1 (mean 0.98).

Moran's I top 10 (NB 34): MGP 0.869, MALAT1 0.763, CRISP3 0.726, CPB1 0.711, IGHG3 0.704,
IGLC2 0.700, IGKC 0.681, ALB 0.668, IGHG1 0.635, IGHG4 0.632. FDR prints as 0.0.
**ERBB2's own Moran's I is never computed or printed.**

Neighbourhood enrichment (NB 36, full 9×9 z matrix): diagonal 51–95; the only positive
off-diagonal pair is 0–1 at z = 8.4; every nest-to-stroma pair runs −4.7 to −22.4.

### 4.6 The per-donor split cannot be re-derived — task stopped

The 3%-to-93% within-tumour range, and HTAPP-313's mean of 2.42, appear only in prose:
NB cell 24 markdown, PDF §2C and §4, DOCX §2.3. No notebook cell computes or prints a per-donor
ERBB2 table, and values for the other six donors are reported nowhere.

**An attempt to re-derive it was stopped at step 1: the AnnData object was never saved.**
The only `.h5ad`/`.h5` files in the module are `data/pbmc3k_raw.h5ad` (the PBMC toy),
`data/anndata/visium_hne_adata.h5ad` (the Visium toy) and
`data/V1_Breast_Cancer_Block_A_Section_1/filtered_feature_bc_matrix.h5` (the spatial matrix, which
carries no cell-type labels or donor IDs). The notebook contains no `write_h5ad` or equivalent
call. Cell 5 pulls the single-cell data live:

```python
census = cellxgene_census.open_soma(census_version="stable")
obs = cellxgene_census.get_obs(census, "Homo sapiens",
        value_filter="disease == 'breast cancer' and is_primary_data == True", ...)
rng = np.random.default_rng(0)
donors = rng.choice(big_donors, size=min(10, len(big_donors)), replace=False)
...
adata = cellxgene_census.get_anndata(census, "Homo sapiens", obs_coords=np.sort(sel), ...)
census.close()
```

Recovering the number would require a fresh network pull of the Census, then re-running QC, because
the per-donor percentages are computed on the 18,388 cells that survive the filters in cell 7
(`filter_genes`, `filter_cells`, `n_genes_by_counts < 6000`, `pct_counts_mt < 15`). Re-running QC
was explicitly out of scope, so the task was not attempted further. Two further obstacles are worth
recording: `census_version="stable"` is a moving pointer, so re-opening it now may resolve to a
release later than the 2025-11-08 one the notebook used, in which case the `soma_joinid` values and
donor sizes differ and the seeded `rng` selects a different donor set; and the ERBB2 values come
from `adata.raw`, since cell 7 subsets the object to 2,000 highly variable genes.

`deliverable_2/lab8_per_donor_erbb2.csv` was therefore not created, and a re-pull has been ruled
out: a moving `stable` Census pointer plus a seeded RNG over a different donor pool produces
different numbers, which is not verification.

**Status of these percentages: `[UNVERIFIED — report prose only]`.** They are not entered in the
provenance log as sourced values, and they may not be cited as measured results. The claim exists in
NB cell 24 markdown, PDF §2C and §4, and DOCX §2.3, and nowhere else. See §7.5 for the write-up
options.

**By contrast, the cell-type contrast that Section 5 rests on IS output-backed.** See §4.4.

### 4.7 COAD content

COAD, colon, CA2 and CR2 appear nowhere in the Lab 8 notebook, PDF, DOCX or README. The transfer
to COAD targets rests entirely on the Lab 7 community result.

---

## 5. Deliverable 1 — Part I targets

Source files: `deliverable_1/{deliverable1.md,evidence_table.csv,inventory.md,provenance.md}`,
`Deliverable1_Bhanushali.pdf`.

Shortlist: CA2 and CR2, with pocket numbers. The full 15-gene triage table carries Lab 2 rank,
Lab 3 rank, week2_score, p_sm, p_ab and final rank. Ranking is lexicographic: integration rank,
then a combined score of 0.6 × omics + 0.4 × prior, then pocket confidence. There is no
druggability threshold.

**The direction problem:** 43 of the top 50 genes, and all 15 candidates, are **down** in tumour
(verified against `targets_coad.csv`). CA2's "honest use is as a differentiation biomarker".
The DOCX handoff names five up-in-tumour top-50 genes to triage next: PLAU, TMEM97, SOX9, DHCR7,
REEP6.

---

## 6. Cross-lab gene map

| Gene | Lab 5 | Lab 6 | Lab 7 | Lab 8 | D1 |
|---|---|---|---|---|---|
| MYC | — | rs6983267, w 0.073 | node, GWAS p 2e-122, PageRank 0.0158 | — | — |
| SMAD7 | — | rs4939827 w 0.111; rs7229639 dropped | node, GWAS p 3e-30 | — | — |
| TCF7L2 | — | rs11196172, w 0.062 | node, GWAS p 7e-20 | — | — |
| BMP4 | — | rs4444235, w 0.086 | node, GWAS p 8e-10 | — | — |
| KRAS | — | — | node, G12D/G12C, PageRank 0.0377 | — | — |
| BRAF | — | — | node, V600E, PageRank 0.0231 | — | — |
| TYMS | DPYD axis (drug, not gene) | — | node, 5-FU + capecitabine, low-conf disease edge | — | — |
| TOP1 | — | — | node, irinotecan, low-conf disease edge | — | — |
| EGFR | — | — | node, PageRank 0.0231 | — | — |
| CA2 | — | — | node, community 9 with Acetazolamide | — | rank 1 |
| CR2 | — | — | node, community 1 (complement) | — | rank 2 |
| DES | — | — | node, community 1 (smooth muscle) | — | top 50 |
| MYH11 | — | — | node, community 1 (smooth muscle) | — | Lab 2 rank 24 |
| DPYD | 105 P/LP, 39 scored | — | — | — | — |
| ERBB2 | — | — | — | 36.4% of tumour cells, mean 0.633 | — |

Genes in two or more Part II labs, for the Phase 2 evidence table: MYC, SMAD7, TCF7L2, BMP4
(Lab 6 + Lab 7); CA2, CR2, DES, MYH11 (Lab 7 + D1). KRAS, BRAF, TYMS, TOP1 and EGFR are Lab 7 only.
DPYD is Lab 5 only. ERBB2 is Lab 8 only. The brief expects rows for all of MYC, SMAD7, TCF7L2,
BMP4, KRAS, BRAF, TYMS, EGFR, CA2 and CR2; six of those will have a single Part II layer and must
be marked as such.

---

## 6a. The tag-versus-causal problem across three labs

The same structure appears in Labs 5, 6 and 7: the variant that was measured is not the variant that
acts, and the distance between them is where the method loses hold. Recorded with the numbers
attached, because this is the connective argument of the case study.

**Instance 1 — Lab 5, DPYD HapB3: the tag is scored, the causal variant is not.**
The exonic tag c.1236G>A (rs56038477) is synonymous, scores **CADD 14.03**, and enters the
calibration as a **benign** variant in band 10–20. The intronic variant that carries the clinical
effect, c.1129-5923C>G (rs75017182), scores **CADD 0.378**, is labelled 3★ drug response, has no
REVEL or AlphaMissense score, and is excluded. The tag is in LD with the variant that causes exon 11
skipping. The calibration sees the tag, takes ClinVar's benign label for it, and never sees the
allele that does the damage.

**Instance 2 — Lab 5, CYP2C19 *2: the variant acts through a mechanism the predictors cannot see.**
rs4244285 is synonymous, so no missense predictor returns a score, and it creates a cryptic splice
donor. It lands at **CADD 5.69** in the 0–10 band, whose LR+ of 0.07 reads as moderate evidence
toward benign. Constraint-based predictors cannot see splicing, and pharmacogene regulatory regions
are under weak constraint by definition, because the selection pressure exists only when the drug
does.

**Instance 3 — Lab 6: 63 tags stand in for unknown causal variants, and the proxy decays with
ancestry.** Each tag was paired to a causal SNP within ±250 kb at r² > 0.6 in EUR. Mean tag-causal
r² then falls from **0.83 in EUR to 0.45 in AFR** (EAS 0.76, SAS 0.77, AMR 0.74;
`module 6/part3_crc/data/README.md`). AFR AUC is 0.567 against EUR 0.685. The score does not
transfer because the tags do not track the same causal variants outside the population whose LD
structure defined them. This is the mechanism behind the equity finding, not a separate fact.

**Instance 4 — Lab 7, rs6983267: the strongest association in the graph sits on a low-confidence
edge, because the gene assignment is a judgement call.** p = **2 × 10⁻¹²²**, the largest in the
graph by 92 orders of magnitude. The `gene` column says MYC; `gwas_mapped_genes` says
CASC8;CCAT2;PCAT1;POU5F1B. A GWAS hit identifies a locus, not the gene it acts through. The edge is
marked low confidence, so the graph is honest about it, but the MYC node exists because of that call.

The four instances differ in what gets lost. Lab 5 loses the mechanism twice over, once to LD
(HapB3) and once to predictor scope (splicing). Lab 6 loses the population, because LD does not
transfer. Lab 7 loses the gene, because annotation was resolved by assumption. In all four the
measured marker stands in for something unmeasured, and in three of them the substitution is silent:
only Lab 7 labels its own uncertainty on the edge.

---

## 7. Decisions needed before Phase 2

### 7.1 The 66 dropped DPYD P/LP variants

DPYD has 105 P/LP variants; the calibration uses the 39 that CADD scored and drops 66 without
comment. The brief's Section 2 plan says "39 of 49 positives are DPYD, so this is effectively a
DPYD calibration". That stands, but the sharper and more defensible statement is that it is a
calibration of the DPYD variants CADD could score, which excludes most frameshift and splice
alleles. Options: state it as a limitation in Section 2, or hold it for Section 7. It is not in
the Lab 5 report either way, so it is new text about existing data, not new analysis.

### 7.2 HapB3 and c.2846A>T

Two clinically actionable DPYD variants are mishandled by the calibration: HapB3's c.1236G>A is
counted as a benign variant in band 10–20, and c.2846A>T is excluded as "Conflicting". Both carry
CPIC activity-score penalties. This is good material for the Section 2 argument about predictor
scope, but it cuts against presenting the DPYD calibration as clean. Include or omit?

### 7.3 The brief's rs7229639 sentence

The planned claim is that "a file defect in one lab quietly weakened a node in another". Verified
facts: rs7229639 is in Lab 6's `score.csv` with a blank effect allele and weight 0.023, and it was
dropped. It is absent from Lab 7. But SMAD7 remains a Lab 7 node with its own PRS tag SNP
(rs4939827, weight 0.111) and a high-confidence GWAS edge at p = 3 × 10⁻³⁰. So the node was not
weakened. Suggested rewrite: the defect removed one of two tags for the SMAD7 locus, and the locus
survived in both labs because a second tag carried it — which makes the point about redundancy
rather than fragility. Confirm the change in framing.

### 7.4 CA2 and the community result — RESOLVED, two separate arguments

**Decision: CA2 is never written into the tissue-composition argument.** Confirmed 2026-11-24.

- **CR2, DES and MYH11** are in community 1 and carry the tissue-composition reading. DES and MYH11
  are smooth-muscle proteins and CR2 is a B-cell and follicular dendritic cell receptor, so the most
  likely interpretation of the gene membership is that the cluster reflects tissue composition
  rather than tumour mechanism: normal colon carries muscle and lymphoid tissue that tumour samples
  largely lack. This is a network-derived hypothesis with no single-cell data behind it.
- **Separately, and not the same weakness:** Constipation and Weight loss also sit in community 1,
  arriving on hand-entered HPO edges marked low confidence, as does the COAD disease node. That is a
  defect in how the cluster was assembled, not evidence about tissue composition. Keep the two points
  distinct in the write-up.
- **CA2** is alone in community 9 with acetazolamide and carbon dioxide hydration. Its weakness is
  **direction**, not cell-type attribution: it is lost in tumour (rna_lfc −5.198, prot_lfc −2.157),
  so inhibiting it reverses nothing. That is Deliverable 1's finding, not a network result.

### 7.5 The Lab 8 per-donor split (affects Section 5)

The 3%-to-93% range has no computed output behind it and cannot be re-derived, because the AnnData
was never saved (§4.6). Under rule 2 it cannot be cited as a measured result. Two options:

1. Drop the range from Section 5 and make the argument from the cell-type table alone, which is
   fully sourced: tumour 36.4% at mean 0.633 against NK 0.258 and everything else at or below 0.06.
   The "a bulk average would blur tumour with stroma" point stands on that table by itself.
2. Report the range explicitly as a claim carried in the Lab 8 report that the notebook does not
   back, and say so in the text.

I recommend 1 for the body and a one-line version of 2 in Section 7, since an unsupported number in
your own prior report is itself a limitation worth naming. Your call.

### 7.5a Druggability values — print convention

Lab 7 (0.769 CA2, 0.640 CR2) and Deliverable 1 (0.770, 0.663) are the **same Module 3 prior, refit
with graph genes held out**, not two measurements. Use Lab 7's values when drawing on Lab 7 and D1's
when citing D1. If both appear in one sentence, say they are the same prior refit.

### 7.6 Repo issues, separate from the deliverable

Not write-up content, but worth fixing in the repo at some point:

- Lab 6 notebook markdown cells 0, 1 and 22 still carry the professor's Alzheimer's template text
  (PGS000334, the APOE ε4 caveat and question) while the analysis and every output are colorectal.
  A grader reading the notebook top to bottom will hit it.
- `modules/module 6/part3_crc/build.log` is from an older build and contradicts the saved data
  (61 of 63 tags, "drop 0", different r² and case rates). It is gitignored, so it does not ship,
  but it is misleading on disk.
- The Lab 7 working-tree notebook is partly un-run and cannot produce the submitted CSVs. The
  committed copy at `bf72daa` is the complete one. Same class of problem as the Lab 5 notebook just
  stashed; worth checking before any future commit.
- The Lab 8 notebook should write the post-QC AnnData to disk (`adata.write_h5ad(...)` after cell 7
  or cell 10) so the object survives the session. As it stands the single-cell results cannot be
  checked, re-plotted or extended without a fresh Census pull that will not reproduce the original
  sample. This is what blocked the per-donor ERBB2 verification in §4.6.

---

## 8. `[MISSING]` register

| Item | Where it was needed |
|---|---|
| GWAS p-values in Lab 6 | Evidence table `GWAS_p` — supplied by Lab 7 `gwas_lookup.csv` instead (§3.6a) |
| Gene annotation for Lab 6 SNPs | Cross-lab mapping — only Lab 7 supplies it |
| rs7229639 in Lab 7 | §7.3 |
| ~~CA2 in the DES/MYH11/CR2 community~~ — resolved: CA2 is in community 9, not a gap but a finding | §7.4 |
| Per-donor ERBB2 table, Lab 8 — AnnData never saved; re-pull ruled out. Marked `[UNVERIFIED — report prose only]`, not logged as sourced | §4.6, §7.5 |
| ERBB2 Moran's I | Figure 4 caption |
| Named single-cell atlas | Figure 4 caption — it is a pooled Census query |
| iLISI, kBET, scIB | Not computed; silhouette only |
| Uncorrected-label ERBB2 table | Backs the "29% of T cells" claim in the Lab 8 report |
| UGT1A1 / irinotecan toxicity data | Section 2 gap — external citation only |
| CPIC DPYD activity-score table | Section 2 — external citation only |
| UC lead SNPs as variant nodes | Only the 4 COAD SNPs are nodes |
| Saved file for the 1.00 Å CR2 side check | D1 disagreement 11 |
| Working-copy outputs, Lab 5 and Lab 7 notebooks | Both read from git HEAD instead |

## 9. Report-versus-data disagreements

Listed, not reconciled. Full list in `sources/`.

**Lab 5.** (1) Graded notebook outputs exist only in git HEAD. (2) Two different "25" sets: cell 7
counts among 80 excluding start-lost, cell 15 among 78 including it; the report treats them as one.
(3) "No missense predictor scored" is true for REVEL, not AlphaMissense (26 blank). (4) Band 25–30
is "moderate, pathogenic, small n" in the notebook but "under 10, not used" in the report.
(5) The *2 G>T row is ignored. (6) "39 DPYD positives" counts only CADD-scored variants; 66 dropped.
(7) HapB3 counted benign, c.2846A>T excluded, neither discussed. (8) The Part 2 notebook's
calibration uses a different filter and different numbers (125 B/LB; LR+ 0.06, 0.08, 2.27, 5.10).

**Lab 6.** (1) `build.log` is from an older build and conflicts with the saved data: 61 of 63 tags,
"drop 0", different r² and case rates. (2) The SAS "about twice" claim holds for Q1 (2.1×) but not
Q2 (1.4×). (3) Notebook markdown cells 0, 1 and 22 still carry the professor's Alzheimer's template
text (PGS000334, APOE ε4) while the analysis is colorectal. (4) The data README says age values are
simulated; the report says "so is the age effect". (5) The README's "one row (none) publishes no
other_allele" is leftover template text. (6) Two report table rows are not notebook outputs: the
overall 10.0% flag rate and the r² row.

**Lab 7 and D1.** (1) The working notebook cannot produce the submitted CSVs. (2) The report calls
CA2, DES, CR2 and MYH11 "my strongest RNA-and-protein genes", but MYH11 is Lab 2 rank 24 and is not
on the 15-gene triage list; ranks 3–10 are skipped silently. (3) Report Table 2 drops JAK2 and IL23R
from the PageRank output. (4) The report's path string substitutes a p-value for the confidence
tokens the notebook printed. (5) Lab 7's `druggability` is the same Module 3 prior refit with graph genes held out
(CA2 0.769118, CR2 0.64031)
while D1 quotes Module 3's (0.770, 0.663) under the same label. (6) Three different CR2–1LY2 RMSDs
circulate: 38.45 Å saved, 2.20 Å in D1, 1.00 Å in `inventory.md` with no saved file. (7) "Eight
confident domains matching CR2's Sushi repeats" versus "eight blocks following CR2's 15 Sushi
repeats". (8) COLEC12 0.131 "next-highest antibody prior" versus "+0.121 next-largest gap".
(9) 3,115 versus 3,099 one paragraph apart. (10) D1 ranks CA2 first while Lab 7 reads its own
cluster as a warning about the list D1 is built on. (11) Two D1 references still flagged `[VERIFY]`.

**Lab 8.** (1) The per-donor split, HTAPP-313's 2.42 and "max count per gene ~250" trace to no
output. (2) "29% of T cells express ERBB2 when uncorrected" has no output behind it. (3) The DOCX
gives 12,121 uncorrected T cells, which is from the 20,000 before QC, not the 18,388 analysed.
(4) Both reports say seven of ten donors came from metastatic-biopsy studies; the legend shows six
HTAPP donors and no output gives study or biopsy site. (5) "CRISP3 marks domains 4 and 7" has no
per-domain CRISP3 numbers. (6) The notebook prints FDR = 0.0; the reports say "FDR < 0.001".
(7) The PDF drops the DCN/LUM/PTPRC ranges and the −7.7 to −9.7 z-scores that the DOCX has, has no
reference list, and contains a typo ("I sawbheld").
