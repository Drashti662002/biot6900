# Multi-omics identification and prioritization of therapeutic targets in colon adenocarcinoma

**BIOT 6900 · Module 2 · Part 3 · Drashti Bharat Bhanushali · 24 September 2026**

*My professor approved swapping Alzheimer's for colon adenocarcinoma. The method is the same one Part 3 asks for: the layers come from different cohorts with no shared patients, so I joined them at the gene level.*

---

## 1. What I built

I put together three layers of real colon cancer data:

| Layer | Cohort | Comparison | n |
|---|---|---|---|
| Transcriptomic | TCGA-COAD (*Nature* 2012;487:330) | tumour vs adjacent normal | 471 vs 41 |
| Proteomic | CPTAC-2 Prospective Colon (Vasaikar, *Cell* 2019;177:1035) | tumour vs adjacent normal, paired | 96 pairs |
| Genomic | The exome data from that same CPTAC study | tumour vs the patient's own blood | 106 patients |

None of these come with a differential table already made, so I calculated each one in `prepare_coad_data.py`. For RNA I used Mann-Whitney, because TCGA tumours and TCGA normals come from different people. For protein I used a Wilcoxon signed-rank test, because CPTAC measured tumour and normal from the same patient and I wanted the extra power that gives. All the URLs, file names and formulas are in `data/PROVENANCE.md`.

Two decisions are worth explaining before any results.

**I checked the protein layer instead of trusting it.** My calculated tumour-minus-normal log2FC matches the log2FC file CPTAC published themselves at Pearson r = 1.0000 across 6,420 shared genes. That only proves I did the arithmetic the same way they did, so I also compared against a separate study — supplementary table S6B of *Cancer Research* 2024;84:2888–2910, a Chinese cohort with 28 paired tumour/normal samples measured label-free rather than by TMT. Across 2,960 shared proteins, 65.5% move the same direction in both cohorts (Pearson r = 0.451), rising to 68.3% if I only look at proteins that were significant in CPTAC. That is a moderate agreement rather than a strong one, which tells me not to lean too hard on any single protein's fold change.

**I had to correct the mutation layer for gene length.** My first attempt used plain mutation frequency and the top of the list was TTN (56.6% of patients), MUC16 (36.8%) and OBSCN (27.4%), all sitting at or above KRAS (33.0%). Those are just long genes collecting passenger mutations. So I estimated one background rate for the whole cohort (1.60 × 10⁻⁵ mutated-gene events per coding base per patient), used each gene's CDS length to work out how often you'd expect it to be hit by chance, and tested the observed count against that with a binomial test. After the correction the top of the genomic layer is TP53 (66.4), APC (50.9), KRAS (43.6), RPL22 (21.5), ACVR2A (17.0), TCF7L2 (14.2), CASP5 (11.9), SOX9 (10.9) and PIK3CA (9.4). That is the expected colorectal driver list, including the genes with coding microsatellites that get hit in MSI tumours. It is a length correction in the spirit of MutSigCV, not a real driver-detection method, and I am not calling it one.

**The join.** All three files use gene symbols. After stripping spaces and upper-casing, the three-way inner join keeps **4,851 genes** out of 17,611 (RNA), 6,709 (protein) and 14,706 (genomic). Proteomics is what costs me the most: RNA and genomics alone would overlap on 10,694 genes. That is why APC, BRAF, MYC, TCF7L2, FBXW7, ACVR2A and MLH1 never appear in my ranking. They were not measured in the proteome, so they could not be scored. Losing APC is the one that bothers me most, since it is mutated in 77.4% of these patients.

**Concordance.** With no shared patients there is nothing to correlate, so concordance is just whether the two layers point the same way. **3,395 of 4,851 genes (70.0%)** agree. If I only keep the 3,115 genes that are significant on both layers, 78.7% agree and 664 genuinely disagree, so most of the disagreement is not just one layer being flat.

---

## 2. Weighting

The usual argument for weighting the genomic layer higher makes sense to me. A somatic mutation sits upstream of everything else, it does not go away when the cell changes state, and stromal contamination cannot create one. RNA and protein levels can be a downstream consequence, or can just reflect a different mix of cell types in the biopsy.

So rather than argue about it, I tested it. I re-scored with weights of 0.25 / 0.25 / 0.50 and measured what happened. The two rankings correlate at Spearman ρ = 0.941 and share 14 of their top 15 genes. KRAS moves from 821 to 445 out of 4,851 and TP53 from 1,392 to 716, but the top of the list barely changes.

I think the reason is that `rank_percentile` throws away most of what makes the genomic layer informative before the weights ever get applied. `neglog10p` runs from 0 to 66 across these genes, but only 69 genes are significant at FDR < 0.05, so once everything is converted to percentiles the vast majority sit in a narrow band and a handful of real drivers get compressed into the same 0–1 range as a gene with one passenger mutation. Doubling the weight on a variable that has already been flattened does not achieve much.

**So I kept equal weights.** Not because they are obviously right, but because I could show that the weights are not the thing holding this analysis back, and changing them would have looked like I had fixed something when I had not. Equal weights also has no free parameters, and I have no independent list of correct answers to tune weights against. Tuning them until familiar genes moved up would just be me fitting to my own expectations.

The things that would actually change this ranking are how each layer's magnitude is defined (|log2FC| versus something based on significance), whether direction is taken into account, and whether percentiles are the right transform for a layer with this much spread.

---

## 3. Top targets

The top 15 under equal weights are CA2, DES, PTN, HSPB6, BCHE, CFL2, GREM2, PRPH, TSPAN7, B3GALT5, CR2, COLEC12, MRGPRF, CP and GNAO1. Every single one is **down** in tumour on both layers.

I do not think that is wrong. Carbonic anhydrase 2 (RNA log2FC −5.20, protein −2.16) and desmin really do drop as tumours lose differentiated colonic epithelium and smooth muscle. But it exposes a problem with scoring on |log2FC|: **43 of my top 50 genes are things the tumour lost.** A gene that disappears from the tumour might be a good diagnostic marker, but you cannot inhibit something that is already gone. Absolute magnitude cannot tell the difference between "the tumour gained something" and "the tumour stopped looking like normal colon".

Reading the same ranking but only looking at genes that went up on both layers is much more useful for picking targets: **PLAU** (rank 22, uPA, involved in invasion), **TMEM97** (27, the sigma-2 receptor, which also has some mutation evidence at neglog10p 3.5), **SOX9** (40), **LRRC15** (56, a fibroblast surface antigen being developed as an antibody-drug-conjugate target), **MMP11** (62), **DPEP1** (85, RNA log2FC +4.95, a membrane-anchored enzyme that has been reported as a colorectal target) and **THBS2** (94).

The known drivers came out at SOX9 (40), CTNNB1 (514), KRAS (821), SMAD4 (895), RPL22 (1,209), TP53 (1,392) and PIK3CA (2,174). **SOX9 is the one gene where all three layers agree**: up at RNA (+2.06), up at protein (+0.58), and mutated in 17.9% of patients (neglog10p 10.9). It is the only classic driver anywhere near the top 50, and it got there by having support on every layer rather than being extreme on one. That is basically the argument for doing multi-omics integration in the first place.

**KRAS came out at rank 821 of 4,851, the 83rd percentile.** I named KRAS at the start as the gene I wanted to investigate, not as a prediction, and I did not change anything to move it. Its profile makes sense once you look at it: KRAS is one of the top three genomic hits (neglog10p 43.6, mutated in 33% of patients) and is almost flat in abundance (RNA −0.44, protein −0.34, significant only because the sample sizes are large). A G12D mutation changes how the protein cycles between GTP and GDP; it does not change how much protein there is. So any score whose RNA and protein components are abundance magnitudes will rank oncogenes that work by point mutation too low, no matter how the layers are weighted. I think that is the single biggest limitation of this scoring scheme, and KRAS happens to show it clearly.

One thing I checked rather than assumed: this is not a G12C cohort. Of 106 patients, 35 (33.0%) have a KRAS coding mutation, and the breakdown is G12D 11, G12V 7, G13D 4, A146T 3, and G12C only 2 (1.9%). The layers I scored are gene-level and carry no allele information at all, so nothing in this analysis is specific to G12C.

---

## 4. A discordant gene: COL5A1

Three genes in my top 60 have `concordant == False`, and two of them are collagens: **COL5A1** (rank 30, RNA **+1.46**, protein **−0.85**, both FDR < 0.05) and **COL1A1** (rank 45, RNA **+2.40**, protein **−1.90**). The third is PACSIN3.

Taking COL5A1: the RNA says collagen V is strongly switched on in tumour, and the protein says there is less of it. The explanation I find most convincing is that the two layers are measuring different things.

- **RNA is a rate, protein is an accumulated amount.** Fibroblasts in the tumour stroma transcribe collagen hard, which is what TCGA sees. But normal bowel wall already has a thick collagen layer built up over years. High synthesis in the tumour does not have to mean more total collagen than in normal tissue.
- **Mature collagen is hard to extract.** It is cross-linked, heavily modified and mostly insoluble, and bottom-up proteomics on tissue lysate systematically under-recovers exactly that fraction. The more cross-linked the matrix, the worse this gets, and tumour matrix is more cross-linked.
- **Collagen gets processed before it is incorporated.** Procollagen has both propeptides cleaved off, so which peptides end up being measured depends on the processing state, and that differs between active tumour stroma and quiet normal matrix.
- **TMT ratios are relative.** Each channel is normalised to total protein, so a tumour with a lot of epithelial cells dilutes every stromal protein including collagens.

What makes me fairly confident this is real rather than one gene being noisy is that two separate fibrillar collagens do the same thing in the same direction, both significant on both layers.

I do have to flag a competing explanation I cannot rule out. Because the RNA and protein layers come from different cohorts, a difference in average tumour purity or stromal content between TCGA and CPTAC could produce this pattern with no post-transcriptional biology involved. Separating the two would need RNA and protein measured on the same tissue, which is what the next section is about.

---

## 5. Limitation: gene-level, cross-cohort

**No patient is in both the RNA cohort and the protein cohort.** So everything here is a statement about genes in colon cancer generally, never about a person. I cannot say that in a given tumour KRAS is mutated and the protein is also up. I can only say those things are separately true of two different groups of people. Any co-occurrence, patient stratification, or claim that two layers agree *within* a tumour is off the table.

I can put a number on what that costs, because CPTAC did run RNA-seq and proteomics on the **same 96 tumours**. Across 7,069 genes the per-gene RNA-protein Spearman correlation has a median of **0.336**, ranging from −0.52 to 0.91; only 24.3% of genes get above 0.5 and 6.4% are negative. Two things follow. First, the Part 1 lesson holds up on real colon data: RNA is not a substitute for protein. Second, my `concordant` column is a very coarse version of that number. "Same direction" compresses a range from −0.52 to +0.91 into True or False, and then reports 70% agreement when the underlying median coupling is only about a third.

That is the contrast with Parts 1 and 2. Matched data lets you correlate per gene within patients, split patients by mutation status, and test post-transcriptional regulation in the same tissue. Unmatched gene-level data gets you a ranked list and nothing finer. The integration level was decided by the data, not by me.

Four other caveats:

1. **Proteomic depth caps everything.** Only 4,851 of roughly 20,000 genes could be scored, and APC never entered the ranking despite being mutated in 77.4% of these patients.
2. **The layers are not comparable designs.** TCGA is unpaired with only 41 normals; CPTAC is 96 matched pairs on TMT ratios. Different normalisation, different power.
3. **My score ignores the p-values I calculated.** Following the method as specified, the magnitudes are |log2FC| and neglog10p, so a big noisy fold change beats a small precise one. I kept the FDR columns in `targets_coad.csv` so Week 3 can filter on them.
4. **The genomic p-values are a simple length correction**, not a validated driver method.

---

## 6. Hand-off

`targets_coad.csv` has all 4,851 genes ranked, with the effect size, p-value and FDR for each layer, the mutation frequency, the `concordant` flag, the equal-weighted `score`, and the genomics-weighted score as a separate column so the sensitivity check is visible downstream.

If I carried this on, the first change would not be the weights. It would be to score direction-aware, so that genes the tumour gained are separated from genes it lost, and to add a codon-level mutation layer instead of a gene-level one, so a gene like KRAS can be ranked on the evidence that actually makes it a target.

To reproduce: `python prepare_coad_data.py`, then Kernel → Restart & Run All in `BIOT6900_Module2_Starter.ipynb`.
