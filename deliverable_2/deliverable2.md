---
title: "Precision medicine in colon adenocarcinoma: what survives when the evidence changes axis"
subtitle: "BIOT 6900 Portfolio Deliverable 2 (Part II Capstone)"
author: "Drashti Bharat Bhanushali"
date: "November 24, 2026"
---

## 1. The disease, the targets, and the question

Colorectal cancer is the third most common cancer worldwide and the second leading cause of cancer death, with about 1.9 million new cases and 904,000 deaths in 2020 (Sung et al., 2021). Colon adenocarcinoma (COAD) is the dominant histology.

Part I nominated two targets, *CA2* and *CR2*, from a 15-gene list built on TCGA-COAD transcriptomics, CPTAC-2 proteomics and a tractability prior. That was somatic evidence: it describes what changed inside a tumour. Part II works mostly on germline variation, which is present from birth, identical in every cell, and describes a person rather than a tumour. A GWAS risk allele and a somatic hotspot are not the same kind of evidence, and this report never treats them as such.

The axis shift reframes the decision as three questions. Who gets screened earlier (Section 3)? Who can safely receive which chemotherapy (Section 2)? Which cells does a target live on (Section 5)? Section 4 sits across all three, because the knowledge graph is the only artifact holding germline, somatic and pharmacologic evidence in one structure. Every number traces to a file or a saved notebook output, listed in `provenance.md`.

## 2. Variant interpretation and pharmacogenomics

### 2.1 *DPYD* changes the dose before the first cycle

*DPYD* encodes dihydropyrimidine dehydrogenase, which clears more than 80% of an administered fluoropyrimidine dose. Carriers of reduced-function alleles accumulate drug and suffer severe, sometimes fatal, toxicity. Fluorouracil and capecitabine are first-line COAD chemotherapy, so this is directly relevant.

The ClinVar set holds 105 pathogenic or likely-pathogenic *DPYD* variants, including the clinically established alleles with their scores: \*2A (rs3918290), splice donor, CADD 23.7, two-star pathogenic, with no REVEL or AlphaMissense score because it is not a missense change; \*13 (rs55886062, I560S), CADD 31.0 (Rentzsch et al., 2019), REVEL 0.92 (Ioannidis et al., 2016), AlphaMissense 0.850 (Cheng et al., 2023); c.2846A>T (rs67376798), CADD 29.9, REVEL 0.962, AlphaMissense 0.293; and HapB3, represented by the synonymous tag c.1236G>A (rs56038477) at CADD 14.03 and the intronic c.1129-5923C>G (rs75017182) at CADD 0.378 (Figure 1C).

CPIC assigns each allele an activity value, sums both alleles into an activity score, and doses from it: avoid fluoropyrimidines at score 0, start at 50% and titrate at score 1 or 1.5, dose normally at score 2 (Amstutz et al., 2018). This carries CPIC evidence level A, and it changes the dose or the drug before the first cycle rather than after a toxic event. Fluorouracil and capecitabine both reach the tumour through thymidylate synthase, appearing as high-confidence ChEMBL edges onto *TYMS* in community 3 of the graph alongside irinotecan and *TOP1*, so the gene governing who tolerates the drug and the target the drug acts on sit in different labs, joined only by the drug name.

### 2.2 The calibration cannot see the variants that matter

I pooled a CADD calibration across five pharmacogenes: 167 labelled variants, 49 pathogenic and 118 benign, binned by CADD with a positive likelihood ratio (LR+, the odds a pathogenic variant falls in a band relative to a benign one) per band (Figure 1B). The top band looks decisive: 34 pathogenic, zero benign, an infinite LR+.

It does not mean what it appears to mean. Of 105 pathogenic *DPYD* variants only 39 carry a CADD score, so 66 are absent from the table: 27 frameshifts, 26 splice variants, 11 with no consequence annotation at all, and 2 others (Figure 1A). CADD scored 27 of 28 stop-gained variants and 7 of 8 missense, but only 5 of 31 splice variants and none of the 27 frameshifts. The failure is representation, not ranking: a scheme built over single-nucleotide substitutions has nothing to assign an insertion or deletion. The eleven unannotated variants, ten duplications and one large deletion, are the strongest version of the point, because the pipeline cannot even say what kind of variant they are, so a reader cannot tell unscoreable from unexamined.

The surviving positive set is 69% stop-gained. The top band separates perfectly because premature stop codons need no model: the consequence annotation already says the protein is truncated. The variants a clinic struggles with, 31 splice and 27 frameshift alleles, are absent from the table that licenses the likelihood ratio. My Module 5 report treated the empty benign cell as a small-sample artifact to correct; the real problem is that it sits in the band holding the easiest positives while the hard ones were excluded by file format. A second exclusion is label noise: c.2846A>T is an established toxicity allele with a CPIC activity value, but ClinVar lists it "conflicting interpretations", so it drops out entirely, and the predictors disagree on it (CADD 29.9 and REVEL 0.962 high, AlphaMissense 0.293). This denominator was chosen twice by the data and never by me, once by file format and once by label.

### 2.3 The scored variant is not the acting variant

HapB3 is the clearest case. Testing reports the synonymous c.1236G>A, CADD 14.03, which ClinVar labels benign, so the calibration counts it as benign in a band whose LR+ of 0.07 reads as evidence toward benign. The allele causing the harm is the intronic c.1129-5923C>G, which creates a cryptic splice site and skips exon 11 (Amstutz et al., 2018), and it scores CADD 0.378. The two are in linkage disequilibrium, so the exonic variant is a tag. The calibration scores the tag, accepts its benign label, and never sees the variant that does the damage. *CYP2C19* \*2 (rs4244285) is the same structure inside one variant: synonymous, so REVEL and AlphaMissense return nothing, acting through a cryptic splice donor, scoring CADD 5.69 in the LR+ 0.07 band, with \*17 (rs12248560) at CADD 3.03 in the same band.

The lesson is scope, not error. Constraint-based predictors cannot see splicing, and pharmacogene regulatory regions are under weak constraint by definition, because the selection pressure only exists when the drug does. Hence an inversion worth stating: *DPYD* is actionable precisely because CPIC rests on functional assays and measured enzyme activity rather than prediction, so the predictors that failed here were never the evidence base for the recommendation. One gap remains: *UGT1A1*, which governs irinotecan toxicity and carries a DPWG recommendation to start poor metabolisers at 70% of the standard dose (Hulshof et al., 2023), was not in my five-gene set, although irinotecan is in the graph through *TOP1*.

Pharmacogenomics answers who can safely receive which drug. The next section asks who should be screened at all.

![**Figure 1. Variant interpretation and the scope of computational predictors.** (A) CADD PHRED scores for the 105 pathogenic or likely-pathogenic *DPYD* ClinVar variants by consequence class. Blue circles are scored; red crosses in the shaded band have no CADD score, shown explicitly rather than omitted. Fractions above each class are scored/total. (B) Positive likelihood ratio by CADD band, pooled across *DPYD*, *CYP2C19*, *CYP2C9*, *TPMT* and *VKORC1* (167 labelled variants). The 25-30 band holds fewer than 10 variants; the 30-50 band holds 34 pathogenic and 0 benign, giving an infinite ratio. (C) Four variants against the bands with ClinVar status. The two HapB3 variants are in linkage disequilibrium: the scored exonic tag is CADD 14.03, the intronic allele that skips exon 11 is CADD 0.378. Data: `modules/module5/data/`.](figures/figure1_variant_interpretation.png){width=90%}

## 3. Polygenic risk: good ranking for one group only

I applied PGS000148, a 63-variant colorectal cancer score from a European-ancestry GWAS (Jeon et al., 2018; Lambert et al., 2021), to 2,336 people from 1000 Genomes phase 3 (Auton et al., 2015). Harmonisation kept 61 variants, flipped 28, and dropped rs7229639 for a blank effect allele.

**The phenotype and the age effect in this cohort are simulated.** Genotypes and weights are real, but case status was generated from a liability model and so was age. Every AUC, calibration curve and observed/expected ratio below measures how a real score behaves on a synthetic outcome. None is evidence about real colorectal cancer risk.

Discrimination, meaning how well the score separates cases from controls, is good in the training ancestry and poor outside it (Figure 2A). The area under the ROC curve is 0.685 in European-ancestry samples and 0.567 in African-ancestry samples, against 0.5 for a coin flip; East Asian 0.689, South Asian 0.606, American 0.670, pooled 0.647.

The mechanism is in the data. Each variant is a tag standing in for an unobserved causal variant nearby, and it only tracks that variant as well as local linkage disequilibrium allows. Mean tag-causal r-squared falls from 0.83 in European-ancestry samples to 0.45 in African-ancestry samples (Figure 2C). African genomes carry shorter haplotype blocks, so a tag chosen in Europeans is a weaker proxy there. The score does not transfer because the tags do not.

Calibration, meaning whether predicted risks match observed rates, tells a more important story (Figure 2B). Fitting a logistic model in the European-ancestry group (intercept -1.054, slope 0.720) and applying it elsewhere, the African-ancestry group has an observed/expected ratio of 0.99, the best average calibration of any non-training group. Yet its curve is nearly flat: predicted risk rises from 0.11 to 0.47 across quintiles while observed risk moves only from 0.16 to 0.34, and the fourth quintile sits below the third. The model gets the group average right while barely ranking anyone within the group, because over-prediction at the top and under-prediction at the bottom cancel.

This is the clearest argument against reporting one summary number: average calibration passed in exactly the group where discrimination failed worst. South Asian shows the reverse, with the average off (O/E 1.20) while ranking stays intact, the bottom quintile under-predicted 2.1-fold (observed 0.32 against predicted 0.15) and the second quintile 1.4-fold.

Flag rates follow (Figure 2C). At a pooled top-10% cutoff of PRS > 1.282, the proportion flagged runs from 6.3% (American) to 14.9% (South Asian). One threshold selects different fractions of different populations, so a screening programme built on it would not deliver equal access. Added to age alone, the score raises cross-validated AUC from 0.579 to 0.666, though both the age effect and the outcome are simulated.

The equity consequence is specific. African-American adults have the highest colorectal cancer incidence and mortality of any racial or ethnic group in the United States (Siegel et al., 2023). The group with the highest real burden is the group this score ranks worst, so deploying it unchanged would widen that gap rather than narrow it.

Risk scores say who to watch. The next section asks which genes the evidence actually supports.

![**Figure 2. Polygenic risk requires discrimination and calibration reported together.** (A) ROC curves and AUCs per ancestry for PGS000148 in 2,336 1000 Genomes participants. (B) Calibration by PRS quintile, from a logistic model fitted in the European-ancestry group and applied to the others; the African-ancestry curve is nearly flat despite an observed/expected ratio of 0.99. (C) Percentage of each group flagged at the pooled top-10% cutoff (PRS > 1.282), with mean tag-causal r-squared in each bar. **Case status and the age effect are simulated**; genotypes and weights are real. Data: `modules/module 6/part3_crc/`.](figures/figure2_polygenic_risk.png){width=100%}

## 4. The knowledge graph: what it supports, and what it only asserts

I built a graph in Neo4j AuraDB (Neo4j, Inc., 2024) of 54 nodes and 61 edges around COAD, drawing genes from my Part I table, drug targets from ChEMBL, pathways from Reactome, phenotypes from HPO and associations from the Catalog lead SNPs of my own PRS, with a source and confidence label on every edge: 25 high, 4 medium, 32 low (Figure 3A).

PageRank ranks *KRAS* first at 0.0377, then *BRAF*, *EGFR* and *TYMS* at 0.0231, ten genes at 0.0158, and *DES*, *MYH11*, *IL23R* and *CR2* at the floor of 0.0086. The caveat belongs next to the ranking: in a graph of this shape a gene gains PageRank almost entirely from inbound drug and variant edges, so the ranking partly measures how many drugs and variants I attached. The evidence is direct: the two-hop query returns 11 approved drug-gene pairs, and *CR2*, *DES* and *MYH11* appear in none of them, having no drug edge at all, which is exactly why they sit at the floor.

Read with that caveat, the ranking does two useful things. First, it recovers established biology unprompted: *KRAS*, *BRAF* and *EGFR* take the top three, and those are the genes COAD treatment actually turns on. Second, it produces one discordant result. *TYMS* ties for second on two drug edges alone, with RNA log fold change +0.85, no protein measurement in my Part I table, and a tractability prior of 0.093. It is a chemotherapy target needed by any dividing cell, not a COAD-specific dependency, and its position reflects the drugs attached to it rather than evidence that it discriminates tumour from normal tissue.

### 4.1 One SNP, two gene assignments

rs6983267 carries the strongest association in the graph, p = 2 x 10^-122^ for colorectal cancer (Fernandez-Rozadilla et al., 2023). The graph attaches it to *MYC* on an edge marked low confidence, and that label is correct, because the GWAS Catalog maps the variant to *CASC8*, *CCAT2*, *PCAT1* and *POU5F1B*, none of which is *MYC* (Figure 3B). A GWAS hit identifies a locus, not the gene it acts through. rs4444235 shows the same pattern: the graph says *BMP4*, the Catalog maps *DDHD1-DT* and *MIR5580*. The cross-lab asymmetry is the finding. My PRS file carries no gene column, so a polygenic score can act on rs6983267 at weight 0.073 without ever committing to a gene. A graph node must be named, so the graph had to commit, and committed to a gene the Catalog does not list, then recorded its own uncertainty on the edge, which is the behaviour to want. *MYC* sits mid-table on PageRank while carrying the strongest p-value in the graph, because centrality never reads the p-value.

### 4.2 The question the graph made easy to ask

Querying for drugs acting on genes of diseases sharing a symptom with COAD returned two candidates, both through ulcerative colitis: ustekinumab via *IL12B* and tofacitinib via *JAK2*. The ustekinumab path runs: ustekinumab targets *IL12B* (ChEMBL, high), *IL12B* associates with ulcerative colitis at p = 1 x 10^-21^ (GWAS Catalog, high), ulcerative colitis has iron deficiency anaemia (HPO, low), COAD has it too (HPO, low). The biology is plausible, since chronic IL-23-driven inflammation is an established route to colitis-associated cancer. The complication is that *IL12B* encodes the p40 subunit shared by IL-12 and IL-23, so ustekinumab blocks both, while the genetic signal implicates IL-23 specifically through *IL23R*; a p19-selective antibody such as mirikizumab would test the hypothesis more cleanly. The weakness is in the edges, not the biology: the two phenotype edges carrying the entire claim were hand-entered by me from HPO and labelled low confidence, and the three shared symptoms are generic to bowel disease. Separately, the *TYMS* and *TOP1* disease edges are low-confidence `pharmacology` edges, present only because a COAD drug targets those genes, making a drug-gene-disease path circular by construction. A graph records assertions with provenance, not truth.

The graph says which genes the evidence supports. The next section asks which cells they live in.

![**Figure 3. The knowledge graph and the gene-assignment problem.** (A) All 54 nodes and 61 edges, redrawn from `nodes.csv` and `edges.csv`. Node colour is entity type, size is PageRank, edge style is confidence (solid high, light medium, dotted low). (B) rs6983267: the strongest association in the graph sits on its weakest kind of edge. The graph asserts *MYC*; the Catalog maps four other genes. Data: `modules/module 7/`.](figures/figure3_knowledge_graph.png){width=100%}

## 5. Cell-type resolution: a method demonstration, then its application

**The single-cell work here is breast cancer, not COAD.** In SCANPY (Wolf et al., 2018) I analysed 18,388 cells from 10 CZ CELLxGENE Census donors (CZI Cell Science Program et al., 2025), integrating with scVI (Lopez et al., 2018) over Harmony (Korsunsky et al., 2019) on silhouette score, plus the 3,798-spot 10x Visium breast cancer section (10x Genomics, 2020) with Squidpy (Palla et al., 2022). The reasoning, not the disease, is what transfers to the COAD targets.

One finding is worth recording for its own sake. The provided labelling function tested for the substring "t cell" before any tumour keyword, and "malignant cell" contains "t cell", so all 8,526 malignant cells were labelled T cells. Had I not checked, every downstream expression statement would have been attributed to the wrong cell type.

With labels corrected, *ERBB2* is confined to tumour cells (Figure 4A): epithelial/tumour cells show mean log-normalised expression 0.633 with 36.4% of 8,678 cells expressing, while every other type sits at or below 0.26 (natural killer cells highest at 0.258). Among the named immune and endothelial populations the fractions expressing are low: T cells 2.5%, macrophages 3.0%, endothelial 3.8%. Dendritic cells (13.9%) and the heterogeneous "Other" group (17.1%) sit higher, but the dendritic population is only 36 cells and "Other" is mixed, so neither changes the picture. Spatially, tumour nests carry mean *ERBB2* of 1.31 to 1.65 against 0.98 and 1.22 in the two stromal and immune domains, and every nest-to-stroma neighbourhood pair sits at z = -4.7 to -22.4, so nests and stroma are segregated rather than interleaved. The Visium donor is independently annotated HER2-positive by 10x Genomics (2020), which is an external check on the expression result rather than a restatement of it.

The general lesson transfers even though the disease does not. A bulk measurement averages tumour, stroma and immune cells, so a gene expressed at 0.633 in a third of tumour cells and near zero elsewhere looks like modest uniform expression. Expression on tumour cells is a different therapeutic proposition from expression on stroma or immune cells: the first is a target, the second a bystander or a liability.

For *CR2*, *DES* and *MYH11*, the graph makes a cell-type claim. Community detection placed all three in one community with the COAD node, smooth muscle contraction and complement cascade. The protein identities point one way: desmin and smooth muscle myosin heavy chain 11 are smooth-muscle structural proteins, and complement receptor 2 sits on B cells and follicular dendritic cells. Normal colon contains a muscularis layer and organised lymphoid tissue, and tumour samples largely replace both, so the most likely reading of all three appearing "lost in tumour" is that the normal samples contained more muscle and lymphoid tissue. That is a statement about what was in the tube, not about what tumour cells are doing. Two cautions belong with it: the reading comes from a network built with no single-cell data, so it is a hypothesis rather than a result; and separately, the same community contains constipation and weight loss on hand-entered low-confidence HPO edges, which is a defect in how the community was assembled rather than evidence about tissue composition.

For *CA2* the problem is different, and it is not a cell-type problem. *CA2* sits alone in its own community with acetazolamide and carbon dioxide hydration, not in the tissue-composition cluster. Its weakness is direction: carbonic anhydrase 2 is down in tumour on both layers, RNA log fold change -5.198 and protein -2.157, so there is no activity to inhibit. An inhibitor removes function the tumour has already lost. That is Deliverable 1's own finding, carried forward unchanged, and it does not depend on cell-type attribution.

The experiment that would settle the *CR2*, *DES* and *MYH11* question is in Section 7.

![**Figure 4. Cell-type and spatial resolution, demonstrated on breast cancer.** (A) *ERBB2* by cell type in 18,388 cells from 10 breast cancer donors (CZ CELLxGENE Census), after correcting a labelling error that had assigned 8,526 malignant cells to the T-cell class. The dashed line marks 0.26, the ceiling for every non-tumour type. (B) Mean *ERBB2* expression against percentage of spots expressing it, by spatial domain, 10x Visium breast cancer section (3,798 spots). Spots hold several cells each, so percentages are uniformly high and the domain means carry the separation. **This dataset is breast cancer, not colon adenocarcinoma**; it demonstrates the method whose reasoning Section 5 applies to the COAD targets. Data: `modules/module 8/`.](figures/figure4_cell_type_resolution.png){width=100%}

## 6. Synthesis: what survives

Table 1 sets out what each Part II layer says about each gene discussed above, including where a layer says nothing.

**Table 1. Layer accounting across the four Part II analyses.** "Not implicated" means the layer was examined and the gene did not appear; "no data" means the layer could not address it.

| Gene | Germline risk (Lab 6) | Pharmacogenomics (Lab 5) | Network (Lab 7) | Cell type (Lab 8 or inferred) | Layers | Status |
|---|---|---|---|---|---|---|
| *DPYD* | Not implicated | 105 P/LP variants; CPIC level A dosing by activity score | Not in graph | No data | 1 + CPIC | **Actionable** |
| *CYP2C19* | Not implicated | \*2 and \*17 both placed in the most-benign band; not a COAD drug gene | Not in graph | No data | 1 | Method finding |
| *MYC* | rs6983267, weight 0.073, p = 2 x 10^-122^ | Not implicated | PageRank 0.0158; gene assignment disputed by Catalog | No data | 2 | Hypothesis |
| *SMAD7*, *TCF7L2* | rs4939827, weight 0.111, p = 3 x 10^-30^; rs11196172, weight 0.062, p = 7 x 10^-20^ | Not implicated | Both PageRank 0.0158; Catalog agrees with both assignments | No data | 2 | Hypothesis |
| *BMP4* | rs4444235, weight 0.086, p = 8 x 10^-10^ | Not implicated | PageRank 0.0158; Catalog maps other genes | No data | 2 | Hypothesis |
| *KRAS*, *BRAF*, *EGFR* | Not implicated | Not implicated | Top three: PageRank 0.0377, 0.0231, 0.0231; G12C/G12D and V600E; sotorasib, adagrasib, encorafenib, cetuximab, panitumumab | No data | 1 | Established, recovered not discovered |
| *TYMS* | Not implicated | Target of the drugs *DPYD* governs | PageRank 0.0231 on two drug edges; disease edge circular | No data | 1 | Discordant |
| *TOP1* | Not implicated | *UGT1A1* absent from my gene set | PageRank 0.0158; irinotecan; disease edge circular | No data | 1 | Gap identified |
| *CA2* | Not implicated | Not implicated | PageRank 0.0158; own community with acetazolamide | Not in the tissue-composition cluster; no single-cell data | 1 | **Not supported: wrong direction** |
| *CR2*, *DES*, *MYH11* | Not implicated | Not implicated | All at PageRank floor 0.0086; no drug edge; one tissue-composition community | Inferred from that community only; no single-cell data | 2 (one inferred) | Hypothesis |

Convergence, and what the absences mean. No gene carries all four layers and most carry one. That is not a failure of the analysis; it follows from the axes. A germline risk locus, a pharmacogene and a somatic drug target are different kinds of object, and asking a polygenic score about *DPYD*, or a pharmacogene panel about *MYC*, is asking the wrong question of the wrong instrument. The absences carry information. *CA2* is absent from the PRS because no common germline variant at that locus associates with colorectal cancer, which tells us its Part I signal was about tumour state rather than inherited risk. *KRAS* and *BRAF* are absent from the Part I list because that list is entirely genes down in tumour, and both are activated rather than lost: the list's construction excluded the two most important genes in COAD treatment.

Who to test. Every COAD patient should receive *DPYD* genotyping before the first fluoropyrimidine dose, covering at minimum \*2A, \*13, c.2846A>T and HapB3. The polygenic score should not be deployed: its outcome here is simulated, it was trained in European-ancestry samples only, and it barely ranks the group with the highest real burden.

This recommendation rests on one Part II layer plus external guideline evidence, and Table 1 makes that visible, so I should address it directly. Pharmacogenomic actionability does not require network corroboration. It requires functional and clinical data: measured enzyme activity, genotype-stratified toxicity rates, prospective trials. CPIC has that evidence base and this analysis does not. A graph edge or a PageRank position would add nothing to the case for *DPYD* testing, and their absence does not weaken it. The single-layer count reflects what my four labs happened to measure, not the strength of the recommendation.

What a genetic risk result means. A PRS flag is a statement about a group average, here on a simulated outcome. It would not entitle a patient to changed treatment and at best informs when screening colonoscopy begins, and its ancestry restriction must be stated to anyone it is applied to, because the score's accuracy depends on their ancestry in a way they cannot be expected to know. A *DPYD* result is categorically different: it changes the dose or the drug, before treatment, for a named mechanism, on a guideline.

Which cells to look at. Cell-type resolution can reverse a target call, because apparent loss may be a change in which cells are present rather than in what tumour cells do. For *CR2*, *DES* and *MYH11* that remains open, and single-cell COAD data would settle it.

What did not survive. *CA2* ranked first in Part I and does not survive. It is lost in tumour on both RNA and protein, so inhibiting it reverses nothing, and the network places it alone with a drug approved for another indication. Saying so is the result: a list built on effect size and tractability, without asking the direction of change or which cells the signal came from, nominated a gene that cannot be drugged in the way the list implied. The polygenic score is not deployable. The IL-23 hypothesis rests on two hand-entered low-confidence edges and three symptoms generic to bowel disease, and the drug that would test it cleanly is not the drug the path returned.

What to do first: order *DPYD* genotyping before fluoropyrimidine chemotherapy, and treat everything else here as a hypothesis awaiting the validation below.

## 7. Limitations and validation

**The Lab 6 phenotype is simulated.** Case status and age were generated from a liability model, so no AUC, calibration curve or O/E ratio in Section 3 is evidence about real colorectal cancer risk.

Ancestry bias is specific. PGS000148 comes from a European-ancestry GWAS. Mean tag-causal r-squared falls from 0.83 to 0.45 between European- and African-ancestry samples, AUC from 0.685 to 0.567, and flag rates at one threshold range from 6.3% to 14.9%.

No claim here describes a patient. Each layer is a different cohort: ClinVar submitters, 1000 Genomes participants, TCGA and CPTAC tumours, ten breast cancer donors. The graph joins them by gene symbol, not by person. Nobody in this analysis has both a *DPYD* genotype and a PRS, so statements combining layers describe genes, not people.

Edges record assertions, not facts. 32 of 61 edges are low confidence. The *TYMS* and *TOP1* disease edges are circular by construction, and the phenotype edges carrying the repurposing path were hand-entered by me.

Predictor scope. CADD, REVEL and AlphaMissense were built on substitutions in constrained coding regions, close to the opposite of where pharmacogenetic variants act. 66 of 105 pathogenic *DPYD* variants have no CADD score, so the calibration describes only the subset that could be scored.

The Lab 8 disease mismatch. The single-cell and spatial work is breast cancer, and its transfer to COAD targets is by reasoning, not measurement. The per-donor *ERBB2* range quoted in my Module 8 report (3% to 93%) could not be verified: no notebook cell computes it, the AnnData object was never saved, and the Census query uses a moving version pointer, so a re-pull would not reproduce the original sample. It is excluded from this report's claims.

The Part I direction problem. 43 of the top 50 Part I genes, and all 15 triage candidates, are down in tumour. An inhibitor cannot restore a lost function, so the list suited biomarker discovery better than target nomination.

Validation plan for the *DPYD* finding. First, confirm the variant calls by an orthogonal method: re-genotype \*2A, \*13, c.2846A>T and HapB3 by Sanger sequencing or a targeted assay on a second sample, excluding array and pipeline artifacts before any dose changes. Second, test the genotype-to-phenotype link: measure dihydropyrimidine dehydrogenase activity in peripheral blood mononuclear cells, or plasma uracil and the uracil/dihydrouracil ratio, and ask whether measured activity matches the CPIC activity score the genotype predicts, which tests the inference rather than the genotype. Third, test clinical benefit prospectively: compare grade 3 or higher toxicity between genotype-guided and standard dosing, with toxicity as the primary endpoint and dose intensity and progression-free survival as secondary endpoints.

For the *CR2*, *DES* and *MYH11* hypothesis: run the Lab 8 pipeline on COAD single-cell data and test whether that expression sits in epithelial cells or in smooth muscle and lymphoid cells. If it sits in muscle and lymphoid cells, the Part I signal was tissue composition and all three should be dropped as targets.

## References

Amstutz, U., Henricks, L. M., Offer, S. M., Barbarino, J., Schellens, J. H. M., Swen, J. J., Klein, T. E., McLeod, H. L., Caudle, K. E., Diasio, R. B., & Schwab, M. (2018). Clinical Pharmacogenetics Implementation Consortium (CPIC) guideline for dihydropyrimidine dehydrogenase genotype and fluoropyrimidine dosing: 2017 update. *Clinical Pharmacology & Therapeutics, 103*(2), 210-216. https://doi.org/10.1002/cpt.911

Auton, A., Brooks, L. D., Durbin, R. M., Garrison, E. P., Kang, H. M., Korbel, J. O., Marchini, J. L., McCarthy, S., McVean, G. A., & Abecasis, G. R. (2015). A global reference for human genetic variation. *Nature, 526*(7571), 68-74. https://doi.org/10.1038/nature15393

Cheng, J., Novati, G., Pan, J., Bycroft, C., Žemgulytė, A., Applebaum, T., Pritzel, A., Wong, L. H., Zielinski, M., Sargeant, T., Schneider, R. G., Senior, A. W., Jumper, J., Hassabis, D., Kohli, P., & Avsec, Ž. (2023). Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science, 381*(6664), eadg7492. https://doi.org/10.1126/science.adg7492

CZI Cell Science Program, Abdulla, S., Aevermann, B., Assis, P., Badajoz, S., Bell, S. M., Bezzi, E., Cakir, B., Chaffer, J., Chambers, S., Cool, J., Covert, M., Gold, B., Goyal, S., Hilton, J., Ji, S., Kiggins, A., Lombardo, M., McCandless, B., ... Weaver, J. (2025). CZ CELLxGENE Discover: A single-cell data platform for scalable exploration, analysis and modeling of aggregated data. *Nucleic Acids Research, 53*(D1), D886-D900. https://doi.org/10.1093/nar/gkae1142

Fernandez-Rozadilla, C., Timofeeva, M., Chen, Z., Law, P., Thomas, M., Schmit, S., Díez-Obrero, V., Hsu, L., Fernandez-Tajes, J., Palles, C., Sherwood, K., Briggs, S., Svinti, V., Donnelly, K., Farrington, S., Blackmur, J., Vaughan-Shaw, P. G., Shu, X.-O., Long, J., ... Dunlop, M. G. (2023). Deciphering colorectal cancer genetics through multi-omic analysis of 100,204 cases and 154,587 controls of European and east Asian ancestries. *Nature Genetics, 55*(1), 89-99. https://doi.org/10.1038/s41588-022-01222-9

Hulshof, E. C., Deenen, M. J., Nijenhuis, M., Soree, B., de Boer-Veger, N. J., Buunk, A. M., Houwink, E. J. F., Risselada, A., Rongen, G. A. P. J. M., van Schaik, R. H. N., Sweep, F. C. G. J., Touw, D. J., van Westrhenen, R., Deneer, V. H. M., Guchelaar, H.-J., & van der Wouden, C. H. (2023). Dutch Pharmacogenetics Working Group (DPWG) guideline for the gene-drug interaction between UGT1A1 and irinotecan. *European Journal of Human Genetics, 31*(9), 982-987. https://doi.org/10.1038/s41431-022-01243-2

Ioannidis, N. M., Rothstein, J. H., Pejaver, V., Middha, S., McDonnell, S. K., Baheti, S., Musolf, A., Li, Q., Holzinger, E., Karyadi, D., Cannon-Albright, L. A., Teerlink, C. C., Stanford, J. L., Isaacs, W. B., Xu, J., Cooney, K. A., Lange, E. M., Schleutker, J., Carpten, J. D., ... Sieh, W. (2016). REVEL: An ensemble method for predicting the pathogenicity of rare missense variants. *The American Journal of Human Genetics, 99*(4), 877-885. https://doi.org/10.1016/j.ajhg.2016.08.016

Jeon, J., Du, M., Schoen, R. E., Hoffmeister, M., Newcomb, P. A., Berndt, S. I., Caan, B., Campbell, P. T., Chan, A. T., Chang-Claude, J., Giles, G. G., Gong, J., Harrison, T. A., Huyghe, J. R., Jacobs, E. J., Li, L., Lin, Y., Le Marchand, L., Potter, J. D., ... Hsu, L. (2018). Determining risk of colorectal cancer and starting age of screening based on lifestyle, environmental, and genetic factors. *Gastroenterology, 154*(8), 2152-2164.e19. https://doi.org/10.1053/j.gastro.2018.02.021

Korsunsky, I., Millard, N., Fan, J., Slowikowski, K., Zhang, F., Wei, K., Baglaenko, Y., Brenner, M., Loh, P., & Raychaudhuri, S. (2019). Fast, sensitive and accurate integration of single-cell data with Harmony. *Nature Methods, 16*(12), 1289-1296. https://doi.org/10.1038/s41592-019-0619-0

Lambert, S. A., Gil, L., Jupp, S., Ritchie, S. C., Xu, Y., Buniello, A., McMahon, A., Abraham, G., Chapman, M., Parkinson, H., Danesh, J., MacArthur, J. A. L., & Inouye, M. (2021). The Polygenic Score Catalog as an open database for the reproducible evaluation and application of polygenic scores. *Nature Genetics, 53*(4), 420-425. https://doi.org/10.1038/s41588-021-00783-5

Lopez, R., Regier, J., Cole, M. B., Jordan, M. I., & Yosef, N. (2018). Deep generative modeling for single-cell transcriptomics. *Nature Methods, 15*(12), 1053-1058. https://doi.org/10.1038/s41592-018-0229-2

Palla, G., Spitzer, H., Klein, M., Fischer, D., Schaar, A. C., Kuemmerle, L. B., Rybakov, S., Ibarra, I. L., Holmberg, O., Virshup, I., Lotfollahi, M., Richter, S., & Theis, F. J. (2022). Squidpy: A scalable framework for spatial omics analysis. *Nature Methods, 19*(2), 171-178. https://doi.org/10.1038/s41592-021-01358-2

Rentzsch, P., Witten, D., Cooper, G. M., Shendure, J., & Kircher, M. (2019). CADD: Predicting the deleteriousness of variants throughout the human genome. *Nucleic Acids Research, 47*(D1), D886-D894. https://doi.org/10.1093/nar/gky1016

Siegel, R. L., Wagle, N. S., Cercek, A., Smith, R. A., & Jemal, A. (2023). Colorectal cancer statistics, 2023. *CA: A Cancer Journal for Clinicians, 73*(3), 233-254. https://doi.org/10.3322/caac.21772

Sung, H., Ferlay, J., Siegel, R. L., Laversanne, M., Soerjomataram, I., Jemal, A., & Bray, F. (2021). Global cancer statistics 2020: GLOBOCAN estimates of incidence and mortality worldwide for 36 cancers in 185 countries. *CA: A Cancer Journal for Clinicians, 71*(3), 209-249. https://doi.org/10.3322/caac.21660

Wolf, F. A., Angerer, P., & Theis, F. J. (2018). SCANPY: Large-scale single-cell gene expression data analysis. *Genome Biology, 19*, 15. https://doi.org/10.1186/s13059-017-1382-0

10x Genomics. (2020, June 23). *Human breast cancer (block A, section 1), Visium spatial gene expression dataset, Space Ranger 1.1.0* [Data set]. https://www.10xgenomics.com/datasets/human-breast-cancer-block-a-section-1-1-standard-1-1-0

Neo4j, Inc. (2024). *Neo4j AuraDB* [Cloud database service]. https://neo4j.com/cloud/aura/
