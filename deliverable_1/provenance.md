# Provenance log: every number in `deliverable1.md`

Paths are relative to the repo root.

**Abbreviations:**
- `T2` = `modules/module2/targets_coad.csv`
- `T3` = `modules/module3/targets_coad_w3.csv`
- `T4` = `modules/module4/structure_assessment_coad.csv`
- `NB2/NB3/NB4` = the executed Lab 2/3/4 notebooks: `modules/module2/BIOT6900_Module2_Starter.ipynb`, `modules/module3/BIOT6900_Module3_Assignment.ipynb`, `modules/module4/BIOT6900_Module4_Assignment.ipynb` (saved cell outputs)
- `ST` = `shared/module4_structures/`
- `ET` = `deliverable_1/evidence_table.csv`, built by `build_evidence_table.py`
- `FIG` = `deliverable_1/make_figures.py`

"Computed" means a trivial count, mean or ratio of an existing column; no new analysis.

## Section 1: Disease and question

| Value | Where | Source |
|---|---|---|
| 1.9 million cases; 904,000 deaths; 2022; 3rd in incidence, 2nd in mortality | ¶1 | Bray et al. (2024), GLOBOCAN 2022 (literature, not repo) |
| 2 of 106 (1.9%) *KRAS* G12C | ¶1 | `modules/module2/data/coad_kras_variant_spectrum.tsv`, row `KRAS_p.G12C`: `n_patients_mutated` = 2, `pct_of_cohort` = 1.89 |
| 11 (10.4%) G12D | ¶1 | same file, row `KRAS_p.G12D`: 11, 10.38 |
| 106 tumours | ¶1 | `T2.n_patients` (106 for every gene) |
| pLDDT ≥ 70 (criterion iii) | ¶1 | Lab 4 rule, NB4 cell 15 (`reject if mean_pLDDT < 70`) |

## Section 2: Multi-omics integration

| Value | Where | Source |
|---|---|---|
| 471 vs 41 samples | ¶1 | `T2.n_tumor` = 471, `T2.n_normal` = 41; `shared/module2_omics/PROVENANCE.md` |
| 96 pairs | ¶1 | `T2.n_pairs` = 96 (max); PROVENANCE.md |
| 106 patients | ¶1 | `T2.n_patients` |
| 4,851 genes | ¶1; Fig 1A title | NB2 cell 25 output "Genes surviving the 3-way join: 4851"; `len(T2)` |
| 17,611 / 6,709 / 14,706 | ¶1 | NB2 cell 25 output (per-layer row counts) |
| 3,950 RNA; 3,776 protein; 26 mutation (FDR < 0.05) | ¶1; Fig 1B | Computed: `(T2.rna_q<0.05).sum()`, `(T2.prot_q<0.05).sum()`, `(T2.gen_q<0.05).sum()` (`FIG` figure1) |
| 3,115 RNA ∩ protein | ¶1 | NB2 cell 27 output "Genes significant on both layers (FDR < 0.05): 3115" |
| 16 all three layers | ¶1; Fig 1B | Computed: `((rna_q<.05)&(prot_q<.05)&(gen_q<.05)).sum()` |
| 3,395 of 4,851 (70.0%) | ¶2 | NB2 cell 27 output "Sign-concordant genes: 3395 / 4851 (70.0%)" |
| 2,451 of 3,115 (78.7%) | ¶2; Fig 1A legend | NB2 cell 27 output "concordant: 2451 (78.7%)" |
| 664 | ¶2; Fig 1A legend | NB2 cell 27 output "discordant: 664" |
| 1,736 not significant on both | Fig 1A legend/caption | Computed: 4,851 − 3,115 (`~both_sig` in `FIG`) |
| 96 CPTAC tumours (matched) | ¶2 | PROVENANCE.md ("across the 96 tumours measured on both") |
| median Spearman 0.336; 7,069 genes | ¶2 | `modules/module2/data/coad_cptac_matched_rna_protein_corr.tsv`: `spearman_r` median = 0.336, rows = 7,069 |
| *COL5A1* +1.46 / −0.85 | ¶2 | `T2`, gene COL5A1: `rna_lfc` 1.4628, `prot_lfc` −0.8465 |
| *COL1A1* +2.40 / −1.90 | ¶2 | `T2`, gene COL1A1: `rna_lfc` 2.4000, `prot_lfc` −1.9031 |
| all 15 top-ranked; mutation q = 0.39–1.00 | ¶3 | `T2.head(15)`: `gen_q` min 0.3913 (CP), max 1.0 (B3GALT5, GNAO1); `ET.n_layers_supporting` = 2 for all 15 |
| *SOX9* rank 40; *KRAS* rank 821 | ¶3 | `T2.rank` for SOX9 (40), KRAS (821) |
| 43 of the top 50 lower in tumour | ¶3 | Computed: `T2.head(50)`, `(rna_lfc<0)&(prot_lfc<0)` = 43 |
| 3,099 (RNA ∩ protein, not mutation) | Fig 1B, caption | Computed: `(R & P & ~G).sum()` in `FIG` |
| 829 / 657 / 0 / 6 / 4 / 240 | Fig 1B | Computed region counts in `FIG` figure1 (printed when the script runs) |
| Top 15 = the 15 rings | Fig 1A | Genes in `T3.gene` |

## Section 3: Druggability triage

| Value | Where | Source |
|---|---|---|
| 15 candidates, all mapped | ¶1 | NB3 cell 4 output "Top 15 targets: 15 mapped, 0 unmapped" |
| ESM C 300M | ¶1 | `shared/module3_pack/manifest.json`: `embedding_models.main.model` = `esmc_300m` |
| AUROC 0.751; AP 0.578; chance 0.333 | ¶1; §6 | NB3 cell 8 output (family split; "Chance level for average precision = 0.333") |
| CA2 *p*~sm~ 0.770 | ¶2; Table 1 | `T3.p_sm` = 0.76984 for CA2 |
| CA1 0.99, CA7 0.96 | ¶2 | `T3.nearest_drug_targets` for CA2 |
| CR2 *p*~ab~ 0.663, *p*~sm~ 0.142 | ¶2; Table 1; §5 | `T3` for CR2: `p_ab` 0.66309, `p_sm` 0.14209 |
| COLEC12 +0.121 (next-largest p_ab − p_sm gap) | ¶2 | `T3.p_ab` − `T3.p_sm` for COLEC12 = 0.13104 − 0.01022 = 0.121. Note: CP has the higher p_ab (0.207) but p_sm > p_ab, so it is not an antibody-over-small-molecule case |
| CSF2RA, IFNGR2, IL1RL2 0.94 | ¶2 | `T3.nearest_drug_targets` for CR2 |
| DES rank 2, PRPH rank 8; 0.001 | ¶3 | `T2.rank` (DES 2, PRPH 8); `T3.p_sm`/`p_ab` 0.00099/0.00137 (DES), 0.00057/0.00075 (PRPH), shown to 3 dp |
| TSPAN7 0.754; BCHE 0.736; CFL2 0.525 | ¶3; §5 | `T3.p_sm` (0.75362, 0.73597, 0.52453) |
| 5.6% antibody-positive training panel | ¶3 | Computed: `shared/module3_pack/panel.tsv`, `ab_clinical` mean over the 1,200-protein training panel (CA2 removed, NB3 cell 6) = 5.6% |

## Section 4: Structural assessment

| Value | Where | Source |
|---|---|---|
| pLDDT < 70 / < 50 thresholds | ¶1 | Lab 4 rule (NB4 cell 15) / AlphaFold convention |
| CA2 mean pLDDT 97.4 | ¶2; Fig 2; Table 1 | `T4.mean_pLDDT`; NB4 cell 7 output; recomputed from `ST/AF-P00918.pdb` CA B-factors = 97.4 |
| 98.5% > 90 | ¶2 | NB4 cell 7 output `pct_above_90` |
| one domain, residues 4–260 | ¶2; Fig 2 | NB4 cell 9 output ("CA2: 1 confident domain(s)", start 4, end 260) |
| 0.37 Å; 257 pairs; 3HS4 | ¶2; §5 | NB4 cell 12 output "257 aligned pairs … RMSD 0.37 Å"; `T4.rmsd_confident` |
| Pocket 6: score 0.317, pLDDT 98.4 | ¶2; Fig 2 | NB4 cell 15 output; `ST/AF-P00918_out/AF-P00918_info.txt` (Pocket 6 Druggability Score 0.317); pLDDT recomputed from `pockets/pocket6_atm.pdb` residues + B-factors = 98.4 (`ET`) |
| Pocket 2: score 0.298, pLDDT 97.9 | ¶2; Fig 2; Table 1; §5 | NB4 cell 15 output; `pocket2_atm.pdb` + B-factors = 97.9 (`ET` notes) |
| 85% of 13 contact residues; pocket 6 covers none | ¶2 | NB4 cell 17 output, AZM A701 row: `site_residues` 13, `site_pocket_covers` 0.85, `best_kept_pocket_covers` 0.00 |
| CR2 mean pLDDT 74.9 | ¶3; Fig 2; Table 1 | `T4.mean_pLDDT`; NB4 cell 7; recomputed from `ST/AF-P20023.pdb` = 74.9 |
| 1.3% > 90; 6.3% < 50 | ¶3 | NB4 cell 7 output (`pct_above_90`, `pct_below_50`) |
| residues 1–19, 1003–1033 | ¶3 | NB4 cell 7 output `disordered_stretches` |
| eight domains; 26.3 Å | ¶3; Fig 2; §5 | NB4 cell 9 output ("8 confident domain(s); mean PAE between domains 1 and 2: 26.3 Å"); `T4.n_domains`, `T4.inter_domain_PAE` |
| Pocket 50: 0.898, pLDDT 65.2 | ¶3; Fig 2; Table 1; §5 | NB4 cell 15 output; `T4` is silent on pocket 50, so from NB4 + `AF-P20023_info.txt` + `pocket50_atm.pdb` B-factors (`ET.top_pocket_score`, `ET.mean_pLDDT_of_pocket_residues`) |
| Sushi 7–8 junction | ¶3; Fig 2 caption | `modules/module4/Module4_Assignment_Report.md` §3 (lining residues 440–515 span Sushi 7–8) |
| Pocket 27: 0.005, pLDDT 75.5 | ¶3; Fig 2; §5 | NB4 cell 15 output; `T4.best_pocket_score` 0.005 |
| 2.20 Å over 80 residues (domain 1) | ¶3 | NB4 cell 12 output "domain 1 (22-147): RMSD 2.20 Å over 80 residues" |
| 38.45 Å | ¶3 | `T4.rmsd_confident` for CR2; NB4 cell 12 output |
| 1LY2 covers residues 20–149 | ¶3 | Module4_Assignment_Report.md §2 (1LY2 DBREF) |

## Section 5: Prioritised target package

| Value | Where | Source |
|---|---|---|
| "2 of 3" layers | Table 1 | `ET.n_layers_supporting` (CA2, CR2 = 2) |
| CA2 −5.20 / −2.16 | Table 1; ¶CA2 | `T2` CA2: `rna_lfc` −5.1977, `prot_lfc` −2.1573 |
| CR2 −1.97 / −1.41 | Table 1; ¶CR2 | `T2` CR2: `rna_lfc` −1.9706, `prot_lfc` −1.4095 (`ET`) |
| CA2 q 2.7 × 10^−22^, 9.5 × 10^−16^ | ¶CA2 | `T2` CA2: `rna_q` 2.725e−22, `prot_q` 9.542e−16 |
| CR2 q 1.7 × 10^−13^, 2.4 × 10^−6^ | ¶CR2 | `T2` CR2: `rna_q` 1.650e−13, `prot_q` 2.413e−06 |
| CA2 5 of 106, q = 0.50 | ¶CA2 | `T2` CA2: `n_mutated` 5, `gen_q` 0.4952 |
| CR2 10 of 106, q = 0.86 | ¶CR2 | `T2` CR2: `n_mutated` 10, `gen_q` 0.8559 |
| integration rank 1 / 11 of 4,851 | ¶ranking; ¶CA2; ¶CR2 | `T2.rank` |
| triage score 1.000 / 0.520; ranks 1 and 8 of 15 | ¶ranking | `T3.final` (1.0, 0.52), `T3.rank` (1, 8) |
| 3–5 targets requested | ¶ranking | The assignment brief |

## Section 6: Limitations and validation

| Value | Where | Source |
|---|---|---|
| 4,851 genes | bullet 3 | as Section 2 |
| APC, BRAF, MYC, TCF7L2, FBXW7, ACVR2A, MLH1 | bullet 3 | NB2 cell 25 output "Known colon genes that didn't make the join" |
| AUROC 0.751 | bullet 5 | as Section 3 |
| 2 of 15 structurally assessed | bullet 7 | `ET.structure_available` = Y for 2 rows |
| 3HS4, pocket 2 | plan step 3 | as Section 4 |

## Added in the 30 October rewrite (decision justifications)

| Value | Where | Source |
|---|---|---|
| Spearman ρ = 0.941; 14 of top 15 shared (0.25/0.25/0.50 weights) | §2 ¶1 | `T2.score` vs `T2.score_genomic_weighted` (Spearman over 4,851 genes; overlap of top-15 by each score); also Module 2 report §2 |
| AUROC 0.833 (random split) | §3 ¶1 | NB3 Step 3 output table ("Random split 0.833 0.707") |
| CA2 removed before fitting (training on 1,200 proteins) | §3 ¶1 | NB3 output "training on 1200 proteins"; `panel.tsv` has 1,201 rows |
| w = 0.6 | §3 ¶1 | Module 3 report §1; NB3 re-ranking cell |
| Reject pLDDT < 70; suspect PAE > 10 Å | §4 ¶1 | Module 4 report §3 rules |
| 15 Sushi repeats in 8 blocks | §4 CR2 | Module 4 report §1 (UniProt Sushi annotation; NB4 domain blocks) |
| 5 up-on-both genes in top 50: PLAU, TMEM97, SOX9, DHCR7, REEP6 | §5 Handoff | `T2` top 50 by `rank` with `rna_lfc` > 0 and `prot_lfc` > 0 (5 genes; COL5A1 and COL1A1 are RNA-up but protein-down) |

## Figures

- **Figure 1:** `FIG` `figure1()` reads `T2` and `T3` only.
- **Figure 2:** `FIG` `figure2()` reads `ST/AF-P00918.pdb`, `ST/AF-P00918_aligned_to_3HS4.pdb` (the Lab 4 alignment, for placing acetazolamide), `ST/3HS4.pdb` (AZM A701 coordinates), `ST/AF-P20023.pdb`, fpocket `pocketN_atm.pdb` files, and `T4` for the mean pLDDT shown in the panel text.

Renders use PyMOL 3.1.4.1 (installed app, not newly installed).
