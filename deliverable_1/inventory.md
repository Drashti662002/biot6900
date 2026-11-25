# Deliverable 1: Phase 1 inventory

Repo root: `~/biot6900`. Disease confirmed as **colon adenocarcinoma (COAD)** in every lab (Module 2 report title, Module 3/4 READMEs, `DISEASE = "coad"` in the Lab 4 notebook). Stop condition not triggered: the Lab 2 ranked table exists.

Paths are relative to the repo root. `modules/*/data` and `modules/module4/structures` are symlinks into `shared/`.

---

## Lab 2: multi-omics integration (`modules/module2/`)

| Item | Location / value |
|---|---|
| Ranked target table (the guide's `targets_ad.csv`) | `modules/module2/targets_coad.csv`: **4,851 rows**, one per gene |
| Columns | `rank, gene, rna_lfc, rna_p, rna_q, mean_tumor_log2tpm, mean_normal_log2tpm, n_tumor, n_normal, prot_lfc, prot_p, prot_q, n_pairs, neglog10p, gen_p, gen_q, mut_freq, n_mutated, cds_len, n_patients, concordant, transcriptomic, proteomic, genomic, score, score_genomic_weighted` |
| Merged multi-omics table | The same file. It *is* the inner-joined table, sorted by `score` |
| Per-layer input tables | `shared/module2_omics/coad_transcriptomics.tsv` (17,611 genes), `coad_proteomics.tsv` (6,709), `coad_genomics.tsv` (14,706) |
| Provenance | `shared/module2_omics/PROVENANCE.md`; built by `modules/module2/prepare_coad_data.py` |
| Notebook | `modules/module2/BIOT6900_Module2_Starter.ipynb` (cells 24–31 = integration, concordance, scoring, export) |
| Report | `modules/module2/Module2_Part3_Report.md` |

**Layers and sources (confirmed from code and PROVENANCE.md):**
- **Transcriptomic:** TCGA-COAD, 471 primary tumours vs 41 adjacent normals; GDC STAR log2(TPM+1) via the UCSC Xena GDC hub; two-sided Mann–Whitney U, BH FDR.
- **Proteomic:** CPTAC-2 Prospective Colon (Vasaikar et al., 2019), 96 paired tumour/normal; PNNL TMT gene-level log ratios via LinkedOmics; Wilcoxon signed-rank, BH FDR.
- **Genomic:** CPTAC-2 colon whole-exome somatic calls (WUSM GATK, gene level) via LinkedOmics, 106 patients. Mutation frequency plus a CDS-length-corrected one-sided binomial test, giving `neglog10p` and `gen_q`. **No CNV data was used.**

**Integration and thresholds actually used:**
- **The only "surviving gene" filter is the three-way inner join on gene symbol: 4,851 genes.** There was no significance or fold-change filter before ranking. Scoring = equal-weight sum of rank percentiles of |RNA log2FC|, |protein log2FC| and genomic −log10 p.
- The only significance threshold anywhere in Lab 2 is **FDR < 0.05**. It was used descriptively for "significant on both layers" (`rna_q`, `prot_q`) and, in the report, for the genomic layer.
- **Concordance (sign agreement, notebook cell 27 output):** 3,395 / 4,851 concordant (70.0%), 1,456 discordant (30.0%). Among 3,115 genes significant on both RNA and protein, 2,451 concordant (78.7%) and 664 discordant. Concordant down: 1,418; concordant up: 1,977.
- **Per-layer FDR < 0.05 counts within the 4,851 joined genes** (computed now from existing `rna_q`/`prot_q`/`gen_q`; a trivial count): RNA 3,950; protein 3,776; genomic 26. Pairwise: RNA∩protein 3,115, RNA∩genomic 22, protein∩genomic 20. **Triple overlap: 16.** None of any layer: 240.
- **Triple-overlap genes:** KCNMA1 (rank 16), COL5A1 (30), DCLK1 (33), SOX9 (40), DDX27 (53), B2M (589), TMA16 (657), NDUFC2 (761), KRAS (821), SMAD4 (895), FXR1 (1037), SEC63 (1040), PHACTR4 (1073), RPL22 (1209), TMEM263 (1706), EGLN1 (1891).

**Integration figure:** none exists. The Lab 2 notebook produces no plots for Part 3, so Figure 1 must be re-plotted from `targets_coad.csv`.

**Top 15 carried to Lab 3:** CA2, DES, PTN, HSPB6, BCHE, CFL2, GREM2, PRPH, TSPAN7, B3GALT5, CR2, COLEC12, MRGPRF, CP, GNAO1. **All 15 are down in tumour on both RNA and protein. All 15 are significant at FDR < 0.05 on RNA and protein, and none is significant on the genomic layer** (`gen_q` 0.39–1.00).

---

## Lab 3: druggability triage (`modules/module3/`)

| Item | Location / value |
|---|---|
| Druggability scores | `modules/module3/targets_coad_w3.csv` (15 rows: `rank, gene, week2_score, p_sm, p_ab, nearest_drug_targets, final`) |
| Embeddings | `shared/module3_pack/proteome_emb_main.npy` (20,190 × 960), `panel_emb_main.npy` (1,201 × 960); baselines `*_emb_baseline.npy` |
| Model | **ESM C 300M** (`esmc_300m`, 960-dim, mean-pooled), per `manifest.json`. Baseline ESM-2 35M (`facebook/esm2_t12_35M_UR50D`); variant scoring ESM-2 650M (`facebook/esm2_t33_650M_UR50D`) |
| Labels | Open Targets tractability, clinical precedence (Approved / Advanced Clinical / Phase 1); UniProt release 2026_03 |
| Classifier | StandardScaler + logistic regression (C = 0.05), one prior per modality: `p_sm` (small molecule), `p_ab` (antibody). The 1 target in the training panel (CA2) was removed before fitting |
| Validation | Family-split (GroupKFold) AUROC 0.751, AP 0.578 (chance 0.333); random split AUROC 0.833, AP 0.707 |
| Re-ranking | `final = 0.6 × pct(Module 2 score) + 0.4 × pct(max(p_sm, p_ab))` |
| Candidates scored | The 15 genes above; all 15 mapped to a reviewed UniProt entry |
| Notebook / report | `BIOT6900_Module3_Assignment.ipynb`; `Module3_Assignment_Report.md` |

**Threshold:** Lab 3 applied **no druggable / not-druggable cut-off**. CA2 and CR2 were chosen for Lab 4 by stated reasoning (CA2 top on both lines of evidence; CR2 the only clear antibody-over-small-molecule split), not by a threshold. See Decision 2 below.

---

## Lab 4: structural assessment (`modules/module4/`)

| Item | Location / value |
|---|---|
| Summary table | `modules/module4/structure_assessment_coad.csv` (CA2, CR2) |
| AlphaFold models | `shared/module4_structures/AF-P00918.pdb` (CA2, 260 aa), `AF-P20023.pdb` (CR2, 1,033 aa); pLDDT in the B-factor column |
| PAE | `AF-P00918-pae.json`, `AF-P20023-pae.json` |
| Experimental structures | `3HS4.pdb` (CA2 + acetazolamide), `1LY2.pdb` (CR2 Sushi 1–2) |
| Pocket detection | **fpocket**; outputs in `AF-P00918_out/`, `AF-P20023_out/` (`*_info.txt` scores, `pockets/pocketN_atm.pdb` lining atoms) |
| Notebook / report | `BIOT6900_Module4_Assignment.ipynb` (Steps 3–7); `Module4_Assignment_Report.md` / `.pdf` |

| | CA2 (P00918) | CR2 (P20023) |
|---|---|---|
| Mean pLDDT | 97.4 | 74.9 |
| % > 90 / % < 50 | 98.5 / 0.4 | 1.3 / 6.3 |
| Confident domains (PAE) | 1 (4–260) | 8; mean PAE between domains 1–2 = 26.3 Å |
| AF vs experimental | 0.37 Å over 257 pairs (3HS4 A) | 38.45 Å in notebook (alignment artefact); 1.00 Å over residues 20–149 in a side check |
| Top fpocket pocket | #6, score 0.317, pocket pLDDT 98.4, keep | #50, score 0.898, pocket pLDDT 65.2, **reject** |
| Drug / ligand site | Acetazolamide (AZM A701) is in **pocket 2** (rank 2, score 0.298, pLDDT 97.9, 85% coverage) | No small-molecule ligand (NAG glycans only) |
| Best kept pocket | #6 (0.317) | #27 (0.005, pLDDT 75.5) |
| Lab 4 call | "Structure supports a small-molecule pocket" | "No confident pocket: consider another modality" (antibody) |

I recomputed the pocket-residue pLDDT values from the saved PDB B-factors and fpocket lining-residue files. All match the Lab 4 notebook to one decimal place.

**Candidates with structures:** CA2 and CR2 only. The other 13 Lab 3 candidates have no Lab 4 structure, pLDDT or pocket assessment.

---

## Gaps and inconsistencies

1. **The genomic layer doesn't support any carried-forward candidate.** All 15 Lab 3 genes have `gen_q` ≥ 0.39. The 16-gene triple overlap contains none of them, so the brief's framing ("the triple-overlap region contains the strongest candidates") is **false for this dataset** and can't be written as stated.
2. **Every candidate is down in tumour on both layers.** The Lab 2 score ranks |log2FC|, so it rewards loss. Labs 2–4 all flag this: "43 of my top 50 are things the tumour lost". It limits any claim that these are inhibition targets.
3. **Only 2 candidates have structures.** A 3–5-target shortlist needs either candidates with no structural evidence or fewer than 3 targets (Decision 3).
4. **No druggability threshold in Lab 3** (Decision 2).
5. **No CNV in the genomic layer.** Mutation only; the `CNV` field is NA throughout.
6. **No Lab 2 integration figure exists**, so it must be re-plotted.
7. **CR2 RMSD of 38.45 Å in `structure_assessment_coad.csv` is an alignment artefact.** The provided global aligner scattered 48 of 130 pairs outside the crystal's range. This is documented in the Lab 4 README and report. The 1.00 Å side-check value was computed outside the notebook and isn't saved to any file; only the report states it.
8. **CA2's top-scoring pocket isn't the drug site.** Pocket 6 ranks first; acetazolamide sits in pocket 2. This must be stated, not glossed over.
9. **fpocket volumes vary slightly between runs** (Lab 4 README). Scores, pLDDT and verdicts don't.
10. **Genomic FDR count differs by denominator.** The Module 2 report's "69 genes significant at FDR < 0.05" counts all 14,706 genomic genes; within the 4,851 joined genes it's 26. Not an error, but the deliverable must use one denominator consistently.
11. **Proteomics coverage removed known COAD drivers from scoring.** APC, BRAF, FBXW7, TCF7L2, ACVR2A, MYC and MLH1 weren't measured in the proteome, so they were never ranked.
12. **ID mapping:** gene symbols throughout Labs 2–3 (TCGA Ensembl IDs mapped to symbols in Lab 2), UniProt accessions from `proteome.tsv` in Lab 3, AlphaFold files named by UniProt in Lab 4. No mismatches for the 15 candidates.
13. **Cross-cohort design.** RNA is TCGA and protein/mutation is CPTAC, with no shared patients, so concordance is gene-level sign agreement only.

---

## Tools available for figures

- **Structure rendering:** PyMOL.app and ChimeraX-1.9.app are installed in `/Applications`, and py3Dmol is in the `biot6900` env. PyMOL can ray-trace PNGs headless, so nothing new needs installing.
- **Venn:** `matplotlib-venn` is not installed. I'd draw the three-set Venn directly with matplotlib (circles plus region counts) rather than installing a package.

---

## Decisions needed before writing

1. **Output format:** DOCX, PDF, or both?
2. **"Druggable (Y/N)":** Lab 3 set no threshold. Either (a) leave it as NA and describe the priors as continuous scores, or (b) adopt an explicit convention for this deliverable, e.g. best-modality prior ≥ 0.5, labelled as a post-hoc reporting cut-off rather than a lab result. Under (b): CA2 0.770, TSPAN7 0.754, BCHE 0.736, CR2 0.663 and CFL2 0.525 are Y; the other 10 are N.
3. **Shortlist size:** only CA2 and CR2 have structural evidence. Options: (a) a **2-target shortlist** (CA2, CR2), stating that the brief asked for 3–5 but structures exist for only 2; or (b) **CA2, CR2 plus 1–3 Lab 3 survivors without structures** (e.g. TSPAN7, BCHE), each flagged "structure not assessed". Rule 1 forbids running new structure or pocket analysis.
4. **Venn thresholds:** I propose FDR < 0.05 on each layer within the 4,851 joined genes, the only threshold used in Lab 2.
5. **The "good target" definition (Section 1):** the brief's definition (multi-layer support, druggable, confident pocket) is met by **no** candidate, because none has genomic support. I'd state that directly and define a good target as: (i) significant and concordant on RNA and protein; (ii) a tractability prior with a stated modality; (iii) a confident structure with a pocket or epitope whose own pLDDT is ≥ 70. The missing genomic layer becomes a named weakness for every candidate, not something hidden.

---

## Decisions (answered 29 October 2026)

1. **Format:** DOCX. Page count is checked through a PDF render.
2. **Druggable (Y/N):** left as NA. The document reports the continuous priors (`p_sm`, `p_ab`) and states that Lab 3 set no threshold.
3. **Shortlist:** **CA2 and CR2 only**, the two candidates with evidence of every type. The document states that the brief asked for 3–5 and explains why only 2 qualify.
4. **Good-target definition:** significant and concordant on RNA and protein; a tractability prior with a stated modality; a pocket or epitope whose own pLDDT is ≥ 70. The missing genomic support is named as a weakness for every candidate.
5. **Venn thresholds:** FDR < 0.05 per layer within the 4,851 joined genes (the only threshold Lab 2 used).
