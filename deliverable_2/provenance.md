# Provenance log — Portfolio Deliverable 2

Every number used in `evidence_table.csv` and destined for the document, with its source file and
column or its saved notebook output. Values marked `[UNVERIFIED]` are claims found in a prior report
that no output supports; they are recorded here so they are not mistaken for sourced values, and
they may not be cited as measured results.

Build script: `build_evidence_table.py`. It reads only the files listed below and writes
`evidence_table.csv`. No analysis is performed: every cell is a lookup, a string join, or a count.

## Primary source files

| Key | Path | What is taken from it |
|---|---|---|
| T | `modules/module 7/data/targets_coad.csv` | 15-gene Part I evidence: `rna_lfc`, `prot_lfc`, `druggability`, `p_sm`, `p_ab`, `rna_q`, `prot_q`, `mut_n`, `mut_q`, `evidence`, `gwas_p`, `gwas_snp`, `gwas_mapped_genes`, `uniprot` |
| GL | `modules/module 7/data/gwas_lookup.csv` | GWAS p-values as written, PMID, study accession, author, Catalog mapped genes |
| E | `modules/module 7/edges.csv` | `src`, `rel`, `dst`, `confidence`, `pvalue`, `source` for all 61 edges |
| N | `modules/module 7/nodes.csv` | 54 nodes with `kind` and node attributes |
| D1 | `deliverable_1/evidence_table.csv` | `lab2_rank`, `druggability_score`, pocket fields, `notes/gaps` for the 15-gene triage list |
| S6 | `modules/module 6/part3_crc/data/score.csv` | `variant`, `effect_allele`, `other_allele`, `weight` (62 rows) |
| R6 | `modules/module 6/part3_crc/data/README.md` | mean tag-causal r² per ancestry; simulated-phenotype statement |
| M5 | `modules/module5/data/DPYD.tsv` (426 rows) and the four other gene TSVs | per-variant `cadd_phred`, `revel`, `alphamissense`, `gnomad_af`, `clinvar_sig`, `clinvar_stars`, `consequence` |
| NB5 | `git show HEAD:"modules/module5/BIOT6900_Module5_Part3_Assignment.ipynb"` | calibration table, coverage counts, VUS counts (working tree had no outputs; see inventory §0) |
| NB6 | `modules/module 6/part3_crc/BIOT6900_Module6_Part3_Graded_CRC.ipynb` | cells 10, 14, 16, 19, 21 outputs |
| NB7 | `git show bf72daa:"modules/module 7/BIOT6900_Module7_Starter.ipynb"` | cells 10, 11, 12, 14, 16 outputs (working tree is partly un-run; see inventory §3) |
| NB8 | `modules/module 8/BIOT6900_Module8_Assignment_Starter.ipynb` | cells 5, 7, 10, 18, 23, 28, 30, 34, 36, 39 outputs |

## evidence_table.csv — column by column

| Column | Source | Notes |
|---|---|---|
| `gene` | selection rule, not a data value | genes in ≥2 of Labs 5–8, or carried from Deliverable 1, plus DPYD and CYP2C19 |
| `axis` | derived from T `evidence` + presence in S6 + presence in M5 | germline = has a PRS SNP; somatic = T `evidence` contains "somatic"; pharmacologic = T `evidence` contains "pharmacology" or the gene is a Lab 5 gene; expression = T `evidence` contains "multi-omics" |
| `part_I_rank` | D1 `lab2_rank` | NA where the gene is not on the 15-gene list. MYH11 carries an explicit note: Lab 2 rank 24, not on the list |
| `in_KG` | N / E | Y if the gene is a node |
| `KG_PageRank` | NB7 cell 14 | 18-row table, 4 dp as printed |
| `KG_community` | NB7 cell 16 | community index plus its membership, verbatim |
| `KG_evidence_type` | T `evidence` + E `TARGETS` rows | drug edges appended as "drug edges: X, Y (ChEMBL, high)" |
| `KG_edge_confidence` | E `confidence` on the gene's `ASSOCIATED_WITH` edge | |
| `in_PRS`, `PRS_SNP`, `PRS_weight` | S6 | weight to 3 dp as published in the file |
| `GWAS_p` | GL `gwas_p`, read as a string | kept verbatim, including the source file's own float artifact on SMAD7 (`3.0000000000000003e-30`) |
| `GWAS_mapped_genes` | GL `gwas_mapped_genes` | populated only where the Catalog mapping excludes the attached gene: MYC and BMP4. NA for SMAD7 and TCF7L2, where the columns agree |
| `PGx_relevance` | M5 + NB5 | DPYD and CYP2C19 only |
| `variant_scores` | M5, per-variant rows | star-allele CADD/REVEL/AlphaMissense values quoted individually |
| `cell_type_evidence` | NB8 cell 23 for ERBB2; NB7 cell 16 for the community inference | states source and disease in every non-NA cell |
| `disease_match` | T `disease`, GL `gwas_trait` | CYP2C19 marked "no": clopidogrel is not a COAD drug |
| `layers_supporting` | derived count over the four Part II layers | germline risk (Lab 6), PGx (Lab 5), network (Lab 7), cell type (Lab 8 or inferred) |
| `provenance_notes` | T, D1, GL | carries rna_lfc/prot_lfc/druggability and the D1-vs-Lab-7 druggability mismatch |

## Values verified directly against source files in this session

| Value | Source | Verified |
|---|---|---|
| DPYD 105 P/LP, 39 with CADD, 66 dropped | M5, recomputed | yes |
| Dropped-set consequences: 27 frameshift, 13 splice acceptor, 13 splice donor, 1 stop, 1 missense | M5, recomputed | yes |
| Scored-set consequences: 27 stop, 7 missense, 3 splice donor, 2 splice acceptor | M5, recomputed | yes |
| *2A rs3918290 CADD 23.7 / 23.4 (two rows) | M5 | yes |
| *13 rs55886062 CADD 31.0, REVEL 0.92, AlphaMissense 0.8504 | M5 | yes |
| c.2846A>T rs67376798 CADD 29.9, REVEL 0.962, AlphaMissense 0.2925, Conflicting | M5 | yes |
| HapB3 tag rs56038477 CADD 14.03, Benign 2★ | M5 | yes |
| HapB3 intronic rs75017182 CADD 0.378, drug response 3★ | M5 | yes |
| Pooled calibration 167 / 49 / 118 and all five band LR+ | M5 recomputed vs NB5 | yes, matches report |
| rs6983267 weight 0.073, effect G, other T | S6 | yes |
| rs4939827 0.111, rs4444235 0.086, rs11196172 0.062 | S6 | yes |
| rs7229639 present with blank `effect_allele`, weight 0.023 | S6 | yes |
| MYC p = 2e-122, PMID 36539618, GCST90129505 | GL row 1 | yes |
| `gene` MYC vs `gwas_mapped_genes` CASC8;CCAT2;PCAT1;POU5F1B, same row | GL row 1 | yes |
| 54 nodes, 61 edges | N, E | yes |
| Fluorouracil→TYMS, Capecitabine→TYMS, Irinotecan→TOP1, Acetazolamide→CA2, all ChEMBL high | E | yes |
| TYMS and TOP1 disease edges are `pharmacology`, confidence low | E | yes |
| CA2 pathway edge is Carbon dioxide hydration only; CA2 in community 9 | E, NB7 cell 16 | yes |
| PageRank: KRAS 0.0377; BRAF/EGFR/TYMS 0.0231; CA2 and 9 others 0.0158; DES/MYH11/IL23R/CR2 0.0086 | NB7 cell 14 | yes |
| All 10 communities and their membership | NB7 cell 16 | yes |
| 11-row 2-hop drug result with status | NB7 cell 10 | yes |
| IL12B path edges: ChEMBL high, GWAS Catalog high, HPO low, HPO low | NB7 cell 12 | yes |
| Per-ancestry n, mean PRS, flag rate, case rate | NB6 cell 14 | yes |
| Per-ancestry AUC | NB6 cell 16 | yes |
| Calibration quintiles, all five ancestries, predicted and observed | NB6 cell 19 | yes |
| EUR-fitted intercept −1.054, slope 0.720 | NB6 cell 19 | yes |
| Age-alone 0.579 vs age+PRS 0.666, 5-fold CV | NB6 cell 21 | yes |
| Mean tag-causal r²: EUR 0.83, AFR 0.45, EAS 0.76, SAS 0.77, AMR 0.74 | R6 | yes |
| ERBB2 by cell type, all 10 types, n / mean / fraction | NB8 cell 23 | yes |
| Integration silhouettes: none 0.095/−0.016, Harmony −0.082/0.006, scVI −0.028/0.035 | NB8 cell 18 | yes |
| 8,526 malignant cells relabelled; Tumor 8,678 / T 2,354 / Other 2,976 | NB8 cell 10 | yes |
| 9 spatial domains; ERBB2 means 1.31–1.65 in nests 2–8 | NB8 cells 28, 30 | yes |
| Moran's I top 10 | NB8 cell 34 | yes |
| Neighbourhood enrichment 9×9 z matrix; only positive off-diagonal 0–1 at z 8.4 | NB8 cell 36 | yes |
| CA2 Part I rank 1, DES 2, CR2 11 | D1 | yes |
| CA2 rna_lfc −5.198, prot_lfc −2.157; CR2 −1.968 / −1.415 | T | yes |

## `[UNVERIFIED]` — claims with no supporting output

| Claim | Where it appears | Status |
|---|---|---|
| ERBB2 detected in 93% of HTAPP-313 tumour cells, 66% of pal_Patient 0176, 3–5% of 0056 and 0178 | NB8 cell 24 markdown; PDF §2C, §4; DOCX §2.3 | `[UNVERIFIED — report prose only]`. No per-donor table is computed anywhere; the AnnData was never saved; a re-pull would not reproduce the sample. Not citable as a measured result |
| HTAPP-313 mean 2.42 | DOCX only | `[UNVERIFIED — report prose only]` |
| "Max count per gene ~250" | DOCX only | `[UNVERIFIED — report prose only]` |
| 29% of T cells express ERBB2 before label correction | PDF, DOCX | `[UNVERIFIED]`. The uncorrected table was never printed |
| Seven of ten donors from metastatic-biopsy studies | PDF, DOCX | `[UNVERIFIED]`. The legend shows six HTAPP donors; no output gives study or biopsy site |
| CRISP3 marks domains 4 and 7 | PDF, DOCX | `[UNVERIFIED]`. No per-domain CRISP3 numbers; the notebook says only "lower-right nests" |
| CR2–1LY2 RMSD 1.00 Å over residues 20–149 | D1 `inventory.md` | `[UNVERIFIED]`. No saved file; 38.45 Å is the saved value and 2.20 Å appears in D1 |

## Values needing an external citation (not in any lab file)

| Value | Needed for | Status |
|---|---|---|
| CPIC DPYD activity-score dosing recommendation and evidence level | Section 2 | external, APA 7, to be added |
| UGT1A1 / irinotecan toxicity and the *28 allele | Section 2 gap paragraph | external; absent from Lab 5 entirely |
| COAD incidence and mortality | Section 1 | external |
| US colorectal cancer burden by ancestry | Section 3 equity paragraph | external |
| HapB3 exon 11 skipping mechanism | Section 2 | external |
| IL-23 p40 vs p19 biology | Section 4 | external; partly discussed in the Lab 7 report |
