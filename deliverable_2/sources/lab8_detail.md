# Lab 8 evidence inventory: ERBB2, single-cell + spatial (breast cancer)

Repo dir: `/Users/drashtibhanushali/biot6900/modules/module 8/`
Sources: `NB` = BIOT6900_Module8_Assignment_Starter.ipynb (cell index 0-based, all 19 code cells executed in order, exec counts 1-19); `PDF` = Module8_ERBB2.pdf (4 pages, no reference list); `DOCX` = Module8_report.docx (longer draft with references); `README` = README.md.
Decoded figure PNGs: `/private/tmp/claude-501/-Users-drashtibhanushali-biot6900-modules-module-10/ee468b5a-f628-4a63-9a5f-4d9be1868906/scratchpad/lab8fig/cell<N>_out<M>.png`

## 1. Datasets

| Item | Value | Provenance |
|---|---|---|
| Single-cell source | CZ CELLxGENE Census, `census_version="stable"` = release 2025-11-08; filter `disease == 'breast cancer' and is_primary_data == True` | NB cell 5 (source + stdout); PDF §1; DOCX §1 |
| Specific atlas/collection name | [MISSING: no named breast-cancer atlas/collection; it is a Census query pooling studies. Donor IDs imply HTAPP, Qian BRCA and "pal_" studies, but no dataset/collection titles are printed] | donor legend in NB cell 10 image |
| Cells matching query | 1,101,430 | NB cell 5 stdout |
| Sampled | 10 whole donors x 2,000 cells = 20,000 cells (rng seed 0; donors with >=2,000 cells) | NB cell 5 |
| After QC | 18,388 cells (filter_genes min_cells=3, min_genes=200, n_genes<6000, pct_mt<15) | NB cell 7 stdout; NB cell 10 stdout |
| Donors | 10: HTAPP-211, HTAPP-313, HTAPP-516, HTAPP-806, HTAPP-889, HTAPP-997, QianBRCA_BC_5, pal_Patient 0056, pal_Patient 0176, pal_Patient 0178 | **image only**: legend of NB cell 10 (cell10_out1.png), also cells 14/18 and PDF Fig 1. Not in any text output |
| Assays (post-QC) | 10x 3' v3 12,635; 10x 3' v2 3,753; 10x 5' v2 2,000 | NB cell 10 stdout; DOCX §1 |
| Visium sample | `V1_Breast_Cancer_Block_A_Section_1` via `sc.datasets.visium_sge` (10x Genomics, Space Ranger 1.1.0 per DOCX refs) | NB cell 26; local copy data/V1_Breast_Cancer_Block_A_Section_1/ |
| Spots x genes | 3,798 x 36,601 | NB cell 26 stdout |
| Spots after `filter_cells(min_counts=500)` | 3,798 (sum of domain n_spots in cell 30 = 3,798; none removed) | NB cell 30 stdout (derived by summing printed n_spots) |
| Other local data | data/pbmc3k_raw.h5ad, data/anndata/visium_hne_adata.h5ad (toy notebooks, not graded) | `ls data` |

## 2. Integration scores (silhouette, 5,000-cell random subsample)

Metric: sklearn `silhouette_score` on embedding, by `donor_id` (batch; lower = better mixing) and by `cell_type_broad` (biology; higher = better). No iLISI/kBET/scIB computed. scVI: n_latent=30, CPU, 60 epochs (not 200). Harmony: harmonypy on 50 PCs, 10 iterations, "Stopped before convergence".

| Method | Batch silhouette | Cell-type silhouette |
|---|---|---|
| No integration (X_pca) | 0.095 | -0.016 |
| Harmony (X_pca_harmony) | -0.082 | 0.006 |
| scVI (X_scVI) | -0.028 | 0.035 |

Provenance: NB cell 18 stdout; NB cell 19 markdown table; PDF Table 1; DOCX Table 1. Chosen method: scVI (NB cell 19, PDF §2B), with the caveat that Harmony is "equally defensible". Harmony log: NB cell 16 stdout. scVI epochs: NB cell 13 stdout ("max_epochs=60 reached").

## 3. Cell-type labels and correction

Broad labels come from keyword function `broad_type()` (NB cell 5). **Correction (NB cell 10):** "malignant cell" contains the substring "t cell", so the T rule matched first. Printed: `malignant cells mislabeled as ['T'] : 8526`. These were reassigned to Epithelial/Tumor before any plot or score (NB cell 10/11; README; PDF §2A; DOCX §2.1).

| Type | Pre-QC, uncorrected (20,000 cells, NB cell 5) | Post-QC, corrected (18,388 cells, NB cell 10) |
|---|---|---|
| Epithelial/Tumor | 152 | 8,678 |
| Other | 3,125 | 2,976 |
| T | 12,121 | 2,354 |
| Endothelial | 1,466 | 1,431 |
| Macrophage | 1,260 | 1,189 |
| Fibroblast | 1,118 | 1,072 |
| B | 586 | 529 |
| Monocyte | 73 | 73 |
| NK | 61 | 50 |
| DC | 38 | 36 |

"Other" = liver/muscle cells from metastasis biopsies (NB cell 19 markdown; DOCX §1, §5). Marker check: CD3D, MS4A1, CD14, EPCAM, PECAM1, COL1A1 on scVI UMAP + dotplot, NB cell 21 (**images only**, no numeric table).

## 4. ERBB2 by cell type (scVI-space labels; log-normalised raw)

NB cell 23 stdout (sorted by mean):

| Type | n_cells | mean_log_expr | frac_expressing |
|---|---|---|---|
| Epithelial/Tumor | 8678 | 0.633 | 0.364 |
| NK | 50 | 0.258 | 0.100 |
| DC | 36 | 0.217 | 0.139 |
| Other | 2976 | 0.192 | 0.171 |
| B | 529 | 0.134 | 0.072 |
| Fibroblast | 1072 | 0.112 | 0.096 |
| T | 2354 | 0.056 | 0.025 |
| Macrophage | 1189 | 0.048 | 0.030 |
| Endothelial | 1431 | 0.042 | 0.038 |
| Monocyte | 73 | 0.000 | 0.000 |

Matches expected: tumour 36% at mean 0.63; all others <=0.26 (max NK 0.258).

**Within-tumour per-donor split** — stated only in prose, NOT produced by any notebook code/output:

| Donor | % tumour cells ERBB2+ | Other | Provenance |
|---|---|---|---|
| HTAPP-313 | 93% | mean 2.42 (DOCX only) | NB cell 24 markdown; PDF §2C; DOCX §2.3 |
| pal_Patient 0176 | 66% | | same |
| pal_Patient 0056 | 3-5% (range given jointly) | | same |
| pal_Patient 0178 | 3-5% (range given jointly) | | same |
| Remaining 6 donors | [MISSING: not reported anywhere] | | |
| Shallow sequencing of low donors | "maximum count per gene about 250" | DOCX §2.3 only | |

[MISSING: the code/output that computed per-donor ERBB2 fractions; no cell groups by donor_id within tumour. Numbers cannot be traced to a saved output.]
Uncorrected-label claim "T cells look like they expressed ERBB2 (29% of them)": PDF §2A, DOCX §2.1. [MISSING: no notebook output shows the uncorrected ERBB2 table.]

## 5. Spatial

Leiden resolution 0.6 on expression only -> **9 domains** (NB cell 28 stdout).

Domain marker means (NB cell 30 stdout, mean log-expression):

| cluster | n_spots | EPCAM | KRT8 | KRT18 | ERBB2 | COL1A1 | COL1A2 | DCN | LUM | PTPRC | CD3D | MS4A1 | CD68 | IGKC |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 725 | 0.66 | 3.72 | 3.36 | 1.22 | 2.34 | 2.37 | 1.19 | 1.20 | 0.44 | 0.17 | 0.08 | 0.66 | 4.20 |
| 1 | 695 | 0.52 | 3.37 | 2.97 | 0.98 | 2.85 | 2.82 | 1.87 | 1.89 | 0.67 | 0.34 | 0.19 | 0.96 | 4.53 |
| 2 | 467 | 0.80 | 4.33 | 3.88 | 1.35 | 2.45 | 2.31 | 0.88 | 1.24 | 0.24 | 0.11 | 0.02 | 0.57 | 2.87 |
| 3 | 412 | 0.93 | 3.74 | 3.33 | 1.38 | 1.78 | 1.67 | 0.71 | 0.87 | 0.24 | 0.12 | 0.04 | 0.53 | 3.01 |
| 4 | 405 | 1.02 | 4.23 | 3.74 | 1.64 | 2.08 | 1.96 | 0.80 | 1.00 | 0.21 | 0.08 | 0.02 | 0.58 | 2.92 |
| 5 | 373 | 1.04 | 3.68 | 3.40 | 1.65 | 1.73 | 1.60 | 0.71 | 0.90 | 0.40 | 0.20 | 0.04 | 0.62 | 2.98 |
| 6 | 358 | 0.99 | 3.47 | 3.20 | 1.31 | 1.53 | 1.51 | 0.56 | 0.67 | 0.15 | 0.05 | 0.03 | 0.42 | 2.93 |
| 7 | 202 | 0.83 | 3.58 | 3.25 | 1.41 | 2.07 | 2.05 | 1.00 | 1.24 | 0.32 | 0.16 | 0.03 | 0.66 | 3.78 |
| 8 | 161 | 1.10 | 3.95 | 3.47 | 1.53 | 1.85 | 1.77 | 0.76 | 0.96 | 0.25 | 0.10 | 0.01 | 0.52 | 2.90 |

Interpretation (NB cell 31; PDF §3A; DOCX §3.1): domains 2-8 = tumour nests; 0-1 = stroma with immune cells.

ERBB2 per domain (NB cell 39 stdout):

| cluster | mean | frac_spots_expressing |
|---|---|---|
| 5 | 1.650 | 0.995 |
| 4 | 1.637 | 0.995 |
| 8 | 1.534 | 1.000 |
| 7 | 1.411 | 0.965 |
| 3 | 1.382 | 0.973 |
| 2 | 1.355 | 0.991 |
| 6 | 1.315 | 0.961 |
| 0 | 1.216 | 0.803 |
| 1 | 0.980 | 0.835 |

**Moran's I** (squidpy `spatial_autocorr`, mode moran; NB cell 34 stdout, top 10 only):

| Gene | I | pval_norm_fdr_bh |
|---|---|---|
| MGP | 0.868504 | 0.0 |
| MALAT1 | 0.762550 | 0.0 |
| CRISP3 | 0.726334 | 0.0 |
| CPB1 | 0.711023 | 0.0 |
| IGHG3 | 0.703817 | 0.0 |
| IGLC2 | 0.700306 | 0.0 |
| IGKC | 0.680696 | 0.0 |
| ALB | 0.667866 | 0.0 |
| IGHG1 | 0.635266 | 0.0 |
| IGHG4 | 0.631525 | 0.0 |

ERBB2 Moran's I: [MISSING: ERBB2 not in top 10 and its I value is never printed or reported in NB, PDF or DOCX].

**Neighbourhood enrichment z-scores** (squidpy nhood_enrichment, 1000 perms; NB cell 36 stdout):

| | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 51.0 | 8.4 | -22.4 | -16.8 | -21.2 | -11.7 | -15.7 | -9.0 | -7.7 |
| 1 | 8.4 | 55.0 | -21.1 | -18.4 | -20.2 | -17.3 | -17.6 | -9.0 | -9.7 |
| 2 | -22.4 | -21.1 | 94.9 | -18.5 | -17.3 | -17.1 | -17.2 | -12.2 | -11.2 |
| 3 | -16.8 | -18.4 | -18.5 | 87.0 | -17.4 | -15.8 | -13.2 | -7.7 | -10.2 |
| 4 | -21.2 | -20.2 | -17.3 | -17.4 | 93.9 | -16.4 | -16.4 | -9.0 | -10.7 |
| 5 | -11.7 | -17.3 | -17.1 | -15.8 | -16.4 | 86.2 | -15.6 | -11.0 | -4.7 |
| 6 | -15.7 | -17.6 | -17.2 | -13.2 | -16.4 | -15.6 | 87.6 | -9.6 | -9.8 |
| 7 | -9.0 | -9.0 | -12.2 | -7.7 | -9.0 | -11.0 | -9.6 | 77.0 | -7.2 |
| 8 | -7.7 | -9.7 | -11.2 | -10.2 | -10.7 | -4.7 | -9.8 | -7.2 | 84.9 |

Readout (NB cell 37; PDF §3C; DOCX §3.3): diagonal 51-95; only off-diagonal positive pair 0-1 (8.4); no tumour-immune interface above chance ("immune-excluded"-like).

## 6. Figure provenance and re-plottability

| Figure | NB cell (decoded PNG) | Report fig | Numbers in text? |
|---|---|---|---|
| Unintegrated UMAP, donor + type | 10 (cell10_out1) | — | No (coords only in image) |
| scVI UMAP donor + type | 14 (cell14_out0) | PDF/DOCX Fig 1A | No |
| Harmony UMAP donor + type | 18 (cell18_out0) | Fig 1B | No; scores yes (cell 18) |
| Labelled map, marker UMAPs, marker dotplot | 21 (out0/1/2) | — | No |
| **ERBB2 on scVI UMAP + cell type** | **23 (cell23_out0)** | **Fig 2** | UMAP itself image-only; **per-type mean / frac / n table is in cell 23 stdout** -> a bar/dot "ERBB2 by cell type" chart can be re-plotted from text without re-running |
| H&E, spatial domains, compartment dotplot | 30 (out1, out2, out3) | — | Domain-mean marker table in cell 30 stdout (dotplot re-plottable as heatmap of means; no fraction-expressing for dotplot sizes) |
| Top-3 Moran genes on tissue (MGP, MALAT1, CRISP3) | 34 (cell34_out2) | — | Moran I table in text; spot-level maps image-only |
| Nhood enrichment heatmap | 36 (cell36_out3) | — | **Full 9x9 z matrix in text** -> re-plottable |
| **Spatial panel: ERBB2 per spot + clusters** | **39 (cell39_out1)** | **Fig 3** | Spot-level maps image-only (needs coordinates); per-domain ERBB2 mean/frac table in cell 39 stdout -> summary chart re-plottable |
| Per-donor ERBB2 within tumour | none | prose only | Not in NB at all; only 4 donors' values in prose |

## 7. COAD / colon / CA2 / CR2

No mention in NB text (sources + text outputs; raw JSON hits for "CA2"/"CR2" are only inside base64 image strings), PDF, DOCX or README. None found.

## 8. Disagreements / untraceable claims (listed, not reconciled)

1. **Per-donor ERBB2 split** (93% HTAPP-313, 66% pal_Patient 0176, 3-5% pal_Patient 0056/0178) appears in NB cell 24 markdown, PDF §2C/§4, DOCX §2.3/§4, but no notebook cell computes or prints it.
2. **HTAPP-313 mean 2.42** and **"maximum count per gene about 250"** for low donors: DOCX only; not in NB or PDF.
3. **"29% of T cells express ERBB2" in uncorrected run** (PDF §2A, DOCX §2.1): no NB output; the uncorrected ERBB2 table is not saved.
4. **"In the uncorrected run this made T the largest group (12,121 cells)"** (DOCX §2.1): 12,121 is the pre-QC 20,000-cell count (NB cell 5), not a post-QC uncorrected count from the 18,388 analysed cells.
5. **"Seven of the ten donors come from metastatic-biopsy (HTAPP) and other atlas studies"** (DOCX §1) / "Seven of the ten donors come from metastatic-biopsy studies" (PDF §1): donor legend (NB cell 10 image) shows 6 HTAPP donors + QianBRCA_BC_5 + 3 pal_Patient; no text output gives donor study/biopsy site. PDF and DOCX also word this differently.
6. **CRISP3 marks "domains 4 and 7"** (PDF §3B, DOCX §3.2) — NB cell 31 only says "lower-right nests"; no per-domain CRISP3 numbers in NB outputs.
7. **Moran's FDR**: NB prints `pval_norm_fdr_bh = 0.0`; reports say "all FDR < 0.001" (consistent in direction, different statement).
8. **PDF vs DOCX content**: the PDF omits DCN/LUM/PTPRC values in §3A, omits z values (-7.7 to -9.7) for the nest 7/8 contacts in §3C, omits the "% of spots" in §3D, and has no reference list (DOCX cites Kedzierska 2025 etc.). PDF has a typo "which then I sawbheld 8,678 cells". PDF has an AI-use statement; DOCX does not.
9. NB cell 23 values agree with reports: 0.63 / 36% tumour; others <=0.26; T 2.5%, macrophage 3%, endothelial 4% (NB 0.038 -> "4%").
10. Integration scores, domain marker ranges, ERBB2 per-domain values and nhood z-scores in PDF/DOCX all match NB outputs.
