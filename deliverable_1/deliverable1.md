---
title: "A target-discovery package for colon adenocarcinoma: from multi-omics evidence to a structurally characterised shortlist"
subtitle: "BIOT 6900 · Portfolio Deliverable 1 (Part I capstone)"
author: "Drashti Bharat Bhanushali"
date: "29 October 2026"
---

# 1. Disease and question

Colorectal cancer accounted for more than 1.9 million new cases and about 904,000 deaths worldwide in 2022, ranking third in incidence and second in cancer mortality (Bray et al., 2024). Targeted options remain confined to molecular subsets. In the CPTAC-2 colon cohort analysed here, only 2 of 106 tumours (1.9%) carry the inhibitor-targetable *KRAS* G12C allele, against 11 (10.4%) with G12D. This package asks which proteins altered in colon adenocarcinoma (COAD) a discovery team could reasonably pursue. A **good target** is defined here by three criteria, each applied in turn in Section 5: (i) significant and direction-concordant change at both RNA and protein level; (ii) a tractability prior that names a plausible modality; and (iii) a confidently predicted structure whose candidate pocket or epitope has its own mean pLDDT ≥ 70. Genomic (somatic mutation) support was intended as a fourth criterion, but no carried-forward candidate meets it, so it is reported as a shared weakness rather than used to rank. The first step is to establish which genes pass criterion (i).

# 2. Multi-omics integration

Three COAD layers were joined at the gene level: TCGA-COAD tumour-versus-normal RNA (471 vs 41 samples; Cancer Genome Atlas Network, 2012; Goldman et al., 2020), CPTAC-2 paired tumour/normal proteomics (96 pairs), and CPTAC-2 exome somatic mutations with a gene-length-corrected recurrence test (106 patients; Vasaikar et al., 2018, 2019). The three-way inner join retained **4,851 genes** (from 17,611 RNA, 6,709 protein and 14,706 mutation-layer genes), all of which were ranked by an equal-weight percentile score. At FDR < 0.05, **3,950** genes are supported by RNA, **3,776** by protein and only **26** by mutation recurrence; **3,115** are supported by both RNA and protein, and **16** by all three layers (Figure 1B).

Direction agreed between RNA and protein for **3,395 of 4,851 genes (70.0%)**, rising to **2,451 of 3,115 (78.7%)** among genes significant on both layers. The remaining **664** genes disagree despite both changes being significant (Figure 1A). Discordance is expected biologically. Within the 96 CPTAC tumours that have both measurements, the median per-gene RNA–protein Spearman correlation is only **0.336** (7,069 genes), reflecting translational control, protein turnover and secretion. The clearest discordant pair is the fibrillar collagens *COL5A1* (RNA log2FC +1.46, protein −0.85) and *COL1A1* (+2.40, −1.90). Their pattern is consistent with stromal fibroblasts transcribing collagen at a high rate while the extractable, cross-linked protein pool remains smaller than in normal bowel wall. Because the RNA and protein layers come from different patients, a difference in tumour purity between cohorts could produce the same pattern, and the two explanations cannot be separated here.

Two features of the ranked list shape everything downstream. First, **all 15 top-ranked genes**, the set carried into triage, are significant and concordant on RNA and protein but **none is significant on the mutation layer** (mutation q = 0.39–1.00). The triple-overlap genes, including *SOX9* (rank 40) and *KRAS* (rank 821), rank lower because the score rewards fold-change magnitude, and oncogenic point mutations such as *KRAS* G12D change protein activity rather than abundance. Second, the score rewards magnitude regardless of sign: **43 of the top 50 genes, and all 15 triage candidates, are lower in tumour than in normal mucosa**. The candidates are therefore strong markers of lost colonic differentiation, and any later claim that one is an *inhibition* target has to overcome that fact. With criterion (i) met by all 15, the next question is which of them a drug could plausibly engage.

![**Figure 1. Multi-omics evidence for 4,851 COAD genes.** **(A)** Tumour-versus-normal log2 fold change at RNA (TCGA-COAD, x-axis) against protein (CPTAC-2, y-axis). Blue circles: direction-concordant genes significant on both layers (FDR < 0.05, n = 2,451). Vermilion triangles: discordant genes significant on both (n = 664). Grey: not significant on both (n = 1,736). Black rings mark the 15 genes carried to druggability triage; all fall in the down/down quadrant. **(B)** Genes significant at FDR < 0.05 in each layer. The 16 genes supported by all three layers (including *SOX9*, *KRAS* and *SMAD4*) contain none of the 15 triage candidates, which all fall in the RNA ∩ protein region (3,099 genes).](figures/figure1_integration_venn.png){width=6.5in}

# 3. Druggability triage

All 15 candidates mapped to reviewed UniProt entries and were scored with two tractability priors, small-molecule (*p*~sm~) and antibody (*p*~ab~). These were logistic-regression classifiers on ESM C 300M sequence embeddings (EvolutionaryScale, 2024), trained against Open Targets clinical-precedence labels (Ochoa et al., 2023). The priors' family-held-out performance (AUROC 0.751, average precision 0.578 against a chance level of 0.333) means they are informative but far from decisive. No druggable/not-druggable threshold was applied; the scores are reported as continuous priors (Table 1; full values in `evidence_table.csv`).

The priors separate the list sharply. **CA2** has the highest small-molecule prior (*p*~sm~ 0.770), and its nearest clinically drugged neighbours in embedding space are CA1 (cosine 0.99) and CA7 (0.96), members of an enzyme family with approved inhibitors (Supuran, 2008). **CR2** is the only candidate whose antibody prior clearly exceeds its small-molecule prior (*p*~ab~ 0.663 vs *p*~sm~ 0.142; the next-highest antibody prior is COLEC12 at 0.131). Its nearest drugged neighbours are the cell-surface cytokine receptors CSF2RA, IFNGR2 and IL1RL2 (all 0.94). These two were carried to structural assessment, as the strongest small-molecule case and the only clear antibody case.

Two results matter for interpreting the shortlist. **DES** (rank 2 in Section 2) and **PRPH** (rank 8) have the lowest priors on the list (*p*~sm~ and *p*~ab~ both 0.001). Both are intermediate-filament proteins with no catalytic site, so strong multi-omics evidence alone does not make a protein druggable. Conversely, **TSPAN7** (*p*~sm~ 0.754), **BCHE** (0.736) and **CFL2** (0.525) have small-molecule priors comparable to or above CR2's best prior, but were not structurally assessed. Their absence from the shortlist reflects unassessed structure, not evidence against them. The antibody prior also rests on few positives (5.6% of the training panel), so CR2's value is the least certain number in this section. Whether either carried-forward protein offers a confident binding site is a structural question.

# 4. Structural assessment

![**Figure 2. AlphaFold models of the two structurally assessed candidates.** Cartoons coloured by AlphaFold pLDDT band (key at bottom); lower panels show per-residue pLDDT with the 70 threshold dashed. **(A)** CA2 (P00918), mean pLDDT 97.4. Purple sticks/points: pocket 2, the acetazolamide site (score 0.298, pocket pLDDT 97.9), with acetazolamide (black) placed from crystal structure 3HS4. Vermilion: pocket 6, the top-scoring pocket (0.317, pLDDT 98.4), which does not overlap the drug site. **(B)** CR2 (P20023), mean pLDDT 74.9. Purple spheres/points: pocket 50, the top-scoring pocket (0.898), whose lining residues average pLDDT 65.2 and dip below 70 at the Sushi 7–8 junction, so it is rejected as a low-confidence pocket claim.](figures/figure2_structures.png){width=6.3in}

AlphaFold models (Jumper et al., 2021; Varadi et al., 2024) were assessed for per-residue confidence (pLDDT), domain-level confidence (predicted aligned error, PAE), agreement with an experimental structure (Berman et al., 2000), and fpocket pockets (Le Guilloux et al., 2009) scored on the mean pLDDT of their own lining residues (Figure 2). Pocket confidence is judged on the pocket region, not the whole protein. A pocket whose residues average pLDDT < 70 is a **low-confidence pocket claim**: side-chain placement and pocket geometry are unreliable, and the region may be disordered or conformationally variable. Below 50, it should not be treated as a pocket at all.

**CA2 (P00918): confident structure, confident pocket.** Mean pLDDT is **97.4**, with 98.5% of residues above 90. PAE resolves a single confident domain (residues 4–260), and the model matches the acetazolamide-bound crystal structure 3HS4 at **0.37 Å** RMSD over 257 residue pairs. The top fpocket pocket (pocket 6, score **0.317**, pocket pLDDT **98.4**) is not the drug site. Acetazolamide sits in **pocket 2** (score **0.298**, pocket pLDDT **97.9**), which covers 85% of the ligand's 13 contact residues. Pocket 6 covers none of them. The pocket that matters is therefore the second-ranked one, and its confidence is as high as the structure's.

**CR2 (P20023): confident domains, uncertain pocket.** Mean pLDDT is **74.9**, with only 1.3% of residues above 90 and 6.3% below 50 (the signal peptide, 1–19, and cytoplasmic tail, 1003–1033). PAE splits the extracellular region into **eight confident domains** matching CR2's Sushi repeats, with a mean inter-domain PAE of **26.3 Å** between the first two, so the relative arrangement of domains is not predicted. The highest-scoring pocket (pocket 50, score **0.898**) has pocket pLDDT **65.2** and lies across the Sushi 7–8 junction. This is a low-confidence pocket claim and plausibly an artefact of two mis-positioned domains. The best pocket above the pLDDT 70 threshold scores only **0.005** (pocket 27, pLDDT 75.5). Sushi 1–2 does agree with crystal structure 1LY2 (2.20 Å over 80 residues for domain 1). The whole-chain RMSD of 38.45 Å reported in the structural table is an alignment artefact, because 1LY2 covers only residues 20–149.

CA2 therefore satisfies criterion (iii) for small molecules, and CR2 does not; CR2 retains a route only through a modality that needs a surface rather than a cavity. The shortlist is built on that distinction.

# 5. Prioritised target package

**Table 1. Ranked target shortlist for COAD.** Layers = number of layers significant at FDR < 0.05 (RNA, protein, mutation). Prior = best-modality tractability prior. Pocket columns refer to the pocket relevant to the modality: the drug site for CA2 and the top-scoring pocket for CR2.

| Rank | Protein | Layers | Concordance (RNA / protein log2FC) | Prior | Mean pLDDT | Pocket pLDDT | Pocket score | Rationale |
|:-------:|:----------|:-------:|:-----------------|:---------|:-------:|:-------:|:---------:|:------------------------|
| 1 | CA2 | 2 of 3 | Concordant, down (−5.20 / −2.16) | *p*~sm~ 0.770 | 97.4 | 97.9 | 0.298 | Confident structure with a known-ligand pocket; small-molecule route |
| 2 | CR2 | 2 of 3 | Concordant, down (−1.97 / −1.41) | *p*~ab~ 0.663 | 74.9 | 65.2 | 0.898 (rejected) | Surface receptor; antibody route; no confident small-molecule pocket |

The shortlist contains two targets rather than the three to five the brief requested. Only CA2 and CR2 carry evidence of all three types (omics, tractability prior and structure); naming further candidates would mean ranking proteins with no structural evidence. Ranking is lexicographic, following the order in which the evidence was generated: integration rank first (CA2 1, CR2 11 of 4,851), then the combined omics-plus-tractability score from triage (CA2 1.000, CR2 0.520; ranks 1 and 8 of 15), then pocket confidence. No tie-break was needed, because CA2 leads on every criterion.

**Rank 1: CA2 (carbonic anhydrase II).** CA2 has the strongest RNA and protein evidence in the dataset: *CA2* RNA log2FC −5.20 (q = 2.7 × 10^−22^) and protein log2FC −2.16 (q = 9.5 × 10^−16^), concordant, giving integration rank 1. Its mutation evidence is weak (*CA2* mutated in 5 of 106 tumours, q = 0.50), the weakness it shares with every candidate. It has the highest small-molecule prior (0.770), its nearest drugged neighbours are other carbonic anhydrases, and it has a structure confident at the pocket itself (pocket pLDDT 97.9, 0.37 Å to a drug-bound crystal). Carbonic anhydrase inhibitors such as acetazolamide are approved drugs (Supuran, 2008), so ligandability is established. CA2 ranks above CR2 because it leads on all three criteria, and because its pocket claim rests on a structure and a ligand position that experiment confirms. **The weakness is biological, not structural:** CA2 is strongly *lost* in tumours, so inhibiting it would not reverse a tumour gain. Its most defensible use is as a differentiation biomarker, or as a test case for whether loss is functionally required for tumour growth.

**Rank 2: CR2 (complement receptor 2, CD21).** *CR2* RNA log2FC −1.97 (q = 1.7 × 10^−13^) and protein log2FC −1.41 (q = 2.4 × 10^−6^) are concordant, giving integration rank 11. Its mutation evidence is also non-significant (10 of 106 tumours, q = 0.86). It is the only candidate with a clear antibody-over-small-molecule prior (0.663 vs 0.142), and its domains are individually well predicted. It is ranked **despite** having no confident small-molecule pocket: the only high-scoring pocket (0.898) sits in a region of pocket pLDDT 65.2 between domains with 26.3 Å PAE, and the best confident pocket scores 0.005. Its route is therefore an antibody against a single well-predicted Sushi domain, not a small molecule. It ranks below CA2 on every criterion. It ranks above BCHE, TSPAN7 and CFL2, which have higher small-molecule priors, only because those three have no structural assessment. Like CA2, CR2 is lower in tumour, which limits it to targeting the cells that still express it.

**Handoff.** We recommend that the discovery team treat CA2 as the lead protein for a validation-first programme (Section 6), with CR2 as a secondary antibody-epitope study, and that structural assessment of TSPAN7 and BCHE precede any expansion of the shortlist.

# 6. Limitations and validation

**What this analysis cannot claim.**

- **Association, not causation.** Lower *CA2* expression in tumours does not show that CA2 loss drives tumour growth.
- **Bulk tissue and cross-cohort data.** Bulk tissue mixes epithelium, stroma and immune cells. The RNA (TCGA) and protein/mutation (CPTAC) layers share no patients, so concordance is gene-level sign agreement and cohort or purity differences cannot be excluded.
- **Proteome coverage.** It limited scoring to 4,851 genes: APC, BRAF, MYC, TCF7L2, FBXW7, ACVR2A and MLH1 were never measured at protein level, so they were never ranked.
- **The mutation layer.** It is a length-corrected recurrence test, not a driver-detection method, and has no copy-number data.
- **Tractability priors.** They encode resemblance to previously drugged proteins (family-split AUROC 0.751), not evidence of binding.
- **Structural limits.** AlphaFold predicts a single static apo conformation, and fpocket on predicted structures yields both false positives (CR2 pocket 50) and mis-ranked true sites (CA2's drug pocket ranks second).
- **Coverage gaps.** Only 2 of 15 candidates were structurally assessed.
- **Direction.** All candidates are decreased in tumour.

**Validation plan for CA2.**

1. **Confirm expression in the right compartment.** Use IHC on matched COAD and normal mucosa, with qPCR and Western blot in COAD cell lines. This establishes whether loss is epithelial or reflects shifts in tissue composition.
2. **Test function in both directions.** Use CRISPR knockout or siRNA knockdown in cell lines that retain CA2, and re-expression in lines that have lost it, with proliferation, viability and differentiation readouts. Cross-check with CRISPR dependency profiles in DepMap (Tsherniak et al., 2017). If CA2 loss is not functionally relevant, the target hypothesis should be dropped in favour of a biomarker use.
3. **Confirm ligandability on the relevant site.** Because pocket confidence is already high and 3HS4 provides an experimental structure, run a thermal-shift (DSF) assay with acetazolamide as a positive control. This confirms that the assay reads out binding in pocket 2 before any new chemical matter is screened.

# References

Berman, H. M., Westbrook, J., Feng, Z., Gilliland, G., Bhat, T. N., Weissig, H., Shindyalov, I. N., & Bourne, P. E. (2000). The Protein Data Bank. *Nucleic Acids Research, 28*(1), 235–242. https://doi.org/10.1093/nar/28.1.235

Bray, F., Laversanne, M., Sung, H., Ferlay, J., Siegel, R. L., Soerjomataram, I., & Jemal, A. (2024). Global cancer statistics 2022: GLOBOCAN estimates of incidence and mortality worldwide for 36 cancers in 185 countries. *CA: A Cancer Journal for Clinicians, 74*(3), 229–263. https://doi.org/10.3322/caac.21834

Cancer Genome Atlas Network. (2012). Comprehensive molecular characterization of human colon and rectal cancer. *Nature, 487*(7407), 330–337. https://doi.org/10.1038/nature11252

EvolutionaryScale. (2024). *ESM Cambrian: Revealing the mysteries of proteins with unsupervised learning* [Technical report]. https://www.evolutionaryscale.ai/blog/esm-cambrian [VERIFY]

Goldman, M. J., Craft, B., Hastie, M., Repečka, K., McDade, F., Kamath, A., Banerjee, A., Luo, Y., Rogers, D., Brooks, A. N., Zhu, J., & Haussler, D. (2020). Visualizing and interpreting cancer genomics data via the Xena platform. *Nature Biotechnology, 38*(6), 675–678. https://doi.org/10.1038/s41587-020-0546-8

Jumper, J., Evans, R., Pritzel, A., Green, T., Figurnov, M., Ronneberger, O., Tunyasuvunakool, K., Bates, R., Žídek, A., Potapenko, A., Bridgland, A., Meyer, C., Kohl, S. A. A., Ballard, A. J., Cowie, A., Romera-Paredes, B., Nikolov, S., Jain, R., Adler, J., … Hassabis, D. (2021). Highly accurate protein structure prediction with AlphaFold. *Nature, 596*(7873), 583–589. https://doi.org/10.1038/s41586-021-03819-2

Le Guilloux, V., Schmidtke, P., & Tuffery, P. (2009). Fpocket: An open source platform for ligand pocket detection. *BMC Bioinformatics, 10*, Article 168. https://doi.org/10.1186/1471-2105-10-168

Ochoa, D., Hercules, A., Carmona, M., Suveges, D., Baker, J., Malangone, C., Lopez, I., Miranda, A., Cruz-Castillo, C., Fumis, L., Bernal-Llinares, M., Tsukanov, K., Cornu, H., Tsirigos, K., Razuvayevskaya, O., Buniello, A., Schwartzentruber, J., Karim, M., Ariano, B., … McDonagh, E. M. (2023). The next-generation Open Targets Platform: Reimagined, redesigned, rebuilt. *Nucleic Acids Research, 51*(D1), D1353–D1359. https://doi.org/10.1093/nar/gkac1046 [VERIFY]

Schrödinger, LLC. (2025). *The PyMOL molecular graphics system* (Version 3.1.4.1) [Computer software]. https://pymol.org

Supuran, C. T. (2008). Carbonic anhydrases: Novel therapeutic applications for inhibitors and activators. *Nature Reviews Drug Discovery, 7*(2), 168–181. https://doi.org/10.1038/nrd2467

Tsherniak, A., Vazquez, F., Montgomery, P. G., Weir, B. A., Kryukov, G., Cowley, G. S., Gill, S., Harrington, W. F., Pantel, S., Krill-Burger, J. M., Meyers, R. M., Ali, L., Goodale, A., Lee, Y., Jiang, G., Hsiao, J., Gerath, W. F. J., Howell, S., Merkel, E., … Hahn, W. C. (2017). Defining a cancer dependency map. *Cell, 170*(3), 564–576.e16. https://doi.org/10.1016/j.cell.2017.06.014

Varadi, M., Bertoni, D., Magana, P., Paramval, U., Pidruchna, I., Radhakrishnan, M., Tsenkov, M., Nair, S., Mirdita, M., Yeo, J., Kovalevskiy, O., Tunyasuvunakool, K., Laydon, A., Žídek, A., Tomlinson, H., Hariharan, D., Abrahamson, J., Green, T., Jumper, J., … Velankar, S. (2024). AlphaFold Protein Structure Database in 2024: Providing structure coverage for over 214 million protein sequences. *Nucleic Acids Research, 52*(D1), D368–D375. https://doi.org/10.1093/nar/gkad1011

Vasaikar, S., Huang, C., Wang, X., Petyuk, V. A., Savage, S. R., Wen, B., Dou, Y., Zhang, Y., Shi, Z., Arshad, O. A., Gritsenko, M. A., Zimmerman, L. J., McDermott, J. E., Clauss, T. R., Moore, R. J., Zhao, R., Monroe, M. E., Wang, Y.-T., Chambers, M. C., … Zhang, B. (2019). Proteogenomic analysis of human colon cancer reveals new therapeutic opportunities. *Cell, 177*(4), 1035–1049.e19. https://doi.org/10.1016/j.cell.2019.03.030

Vasaikar, S. V., Straub, P., Wang, J., & Zhang, B. (2018). LinkedOmics: Analyzing multi-omics data within and across 32 cancer types. *Nucleic Acids Research, 46*(D1), D956–D963. https://doi.org/10.1093/nar/gkx1090
