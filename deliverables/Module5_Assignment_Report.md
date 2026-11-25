# Variant interpretation in CYP2C19

**BIOT 6900 · Module 5 Assignment · Drashti Bharat Bhanushali · 23 October 2026**

*My gene is CYP2C19, the liver enzyme that turns clopidogrel into its active form and breaks down drugs like voriconazole, PPIs and some SSRIs. The fetch gave me 80 ClinVar variants. Both anchor alleles were there: \*2 (rs4244285, two rows, G>A and G>T at the same position) and \*17 (rs12248560). CADD scored 85% of them, REVEL and AlphaMissense 67.5% each, and 68.8% had a gnomAD frequency. I did not have an AlphaGenome key, so I skipped the Atlas (Step 7) and answered that part from reasoning.*

**The one-line version of what I found:** the tools in this course are built to answer "does this variant cause disease", and CYP2C19 does not cause disease. Every problem below is a version of that mismatch.

---

## 1. Thresholds and calibration

**Frequency threshold: 0.5.** In a rare-disease gene you can throw out anything above 0.1%, because a variant that causes a serious disease can't be common. That logic doesn't work for CYP2C19. Losing CYP2C19 function doesn't make anyone sick. It only matters when the person takes one of its drugs, so there has been nothing to keep these alleles rare. \*2 is at 16.8% in my data and \*17 at 20.0%. So the filter has to keep common variants. At 0.5 the only things removed are variants where the "alt" allele is actually what most people carry, meaning the reference genome has the rarer one. It removed 2 of 80 variants, both benign: rs17885098 (synonymous, AF 0.919) and rs3758581 (I331V, AF 0.952). No P/LP, drug-response or anchor variant was lost. I kept the 25 variants that have no gnomAD frequency, because "not in gnomAD" doesn't mean "common".

**The 0.001 contrast.** A rare-disease filter would have removed 19 variants, and 14 of them are drug-response or anchor variants. That list includes \*2 (rs4244285), \*17 (rs12248560), the stop-gain W212X (rs4986893) and the start-loss M1V (rs28399504). In other words, it would have deleted the alleles that CYP2C19 testing is actually done for.

**Predictor: CADD.** It is the only tool that scores almost every class in my gene (68/80). REVEL and AlphaMissense only see missense changes, and \*2 isn't one.

**CYP2C19 could not be calibrated on its own.** After filtering it had **0 P/LP** and 5 B/LB variants with a CADD score. ClinVar files CYP2C19 alleles as "drug response", not pathogenic, so there is nothing to count as a positive. I pooled it with the other four genes from the lab (CYP2C9, VKORC1, DPYD, TPMT), using the same fetch and the same 0.5 filter.

| CADD band | n | P/LP | B/LB | LR+ | Evidence |
|---|---|---|---|---|---|
| 0–10 | 76 | 2 | 74 | 0.07 | moderate, benign |
| 10–20 | 34 | 1 | 33 | 0.07 | moderate, benign |
| 20–25 | 17 | 8 | 9 | 2.14 | supporting, pathogenic |
| 25–30 | **6** | 4 | 2 | 4.82 | **under 10, not used** |
| 30–50 | 34 | 34 | 0 | inf | strong, pathogenic |

The whole calibration rests on **167 labelled variants (49 P/LP, 118 B/LB)**. That doesn't make me very confident, for four reasons:

- **It's really a DPYD calibration.** 39 of the 49 positives come from DPYD, and CYP2C19 adds only 5 labelled variants, all benign.
- **One of those 5 benign labels is \*17**, which ClinVar marks "Benign; other". \*17 is an increased-function allele, so even the benign side for my gene isn't clean.
- **Two bands have problems.** The 25–30 band rests on 6 variants, so I didn't apply it. The 30–50 band has zero benign variants, so its LR is infinite. Adding 0.5 to every cell gives about 165, which is still strong, but that number comes from an empty cell, not a measurement.
- **The benign bands would misclassify \*2.** The 0–10 band says "moderate, benign" (LR 0.07), and \*2 has a CADD score of 5.7, so it falls right in that band. \*17 (CADD 3.0) does too. A calibration that would call my gene's most important no-function allele "moderate benign" can't be trusted for CYP2C19 on the benign side.

**VUS results.** Of 33 VUS, 8 lean pathogenic (all in the 20–25 band, supporting strength only), 16 lean benign, and 9 stay put: 3 sit in the 25–30 band, and 6 are large deletions with no score. Given the problems above, I treat the 16 "leaning benign" calls as the weakest result in this report. None of the 33 calls is P/LP or B/LB, because one computational code alone can't reclassify a variant.

---

## 2. Coverage gaps and the regulatory genome

**Classes no missense predictor scored.** REVEL returned nothing for 8 classes, covering 25 of my 80 variants:

| Class | n |
|---|---|
| Large deletion / CNV (no consequence annotation) | 8 |
| Synonymous | 7 |
| Upstream (promoter) | 5 |
| Intergenic | 1 |
| Intronic | 1 |
| Splice donor | 1 |
| Stop gained | 1 |
| Stop lost | 1 |

The two that bother me most are the **stop_gained** and the **splice_donor**. Those are severe by any standard, and REVEL and AlphaMissense are blank on them only because they aren't amino-acid substitutions. A pipeline reading REVEL alone would show the same empty cell for a nonsense variant as for a variant it had never heard of, and nothing in the output would say which was which.

**My unscored variant: rs72558186 (CYP2C19\*7), chr10:g.94781999T>A.** This is the clearest case in my gene. It is a splice donor variant at the exon 5 / intron 5 boundary, ClinVar lists it as drug response with a **4-star (practice guideline)** review status, and every column I have is empty: no CADD, no REVEL, no AlphaMissense, and no gnomAD frequency. It tripped `flag_no_tool_at_all`. So the single variant in my gene with the strongest clinical review status is the one my pipeline knows least about.

**What evidence would classify it, and what would produce that evidence.** Because it sits in the donor motif, the claim to test is about RNA, not protein, so no missense tool can ever contribute:

1. **SpliceAI**, for the donor-loss delta score. This is the fast first step and it is the real version of the splice model I built in Part 1. A donor-loss delta above 0.8 is the high-precision threshold.
2. **A minigene splicing assay, or RNA-seq from carrier liver tissue**, for the definitive answer. Clone the exon with its flanking intron, transfect, then RT-PCR the product and look for exon skipping or intron retention. SpliceAI predicts the effect; this observes it.

**The Atlas, which I skipped.** I had no API key, so `SKIP_ATLAS = True`. The cost is specific rather than general: the 25 non-missense variants above are exactly what the Atlas is for, and for three of them it is the only tool in this course that could have said anything.

- **rs12248560 (\*17), AF 0.200, CADD 3.03.** A promoter variant that increases transcription and is used clinically. CADD's 3.03 is a shrug, because a promoter under weak evolutionary constraint doesn't look remarkable on a genome-wide deleteriousness scale. The Atlas could have returned a `MAX_ABS_DNASE` or `CHIP_TF` attribution naming a tissue and a chromatin effect, which is a mechanism rather than a constraint score.
- **rs12777823, AF 0.190, CADD 0.725.** Intergenic, drug response, and essentially invisible to everything I ran.
- **rs72558186 (\*7).** A `MERGED_SPLICING` attribution would have flagged the splicing mechanism independently of SpliceAI.

**Why nothing else in the course can do this.** CADD, REVEL and AlphaMissense all ultimately read evolutionary constraint, and pharmacogene regulatory regions are under weak constraint by definition, because the selection pressure only exists when the drug does. AlphaGenome predicts molecular readouts from sequence, per tissue, which is a different signal class. That also means it is the only tool here that could give me a genuinely independent second opinion, since the Part 2 correlation check showed my three constraint-based scores are close to one opinion stated three times.

**The caveat I would still have applied.** A high AVI driven by `CACTUS_241_WAY` or `PHASTCONS_470_WAY` means conservation is carrying the score and the mechanism is unresolved, which for a pharmacogene is close to no information. I would have read the attribution column before the number.

---

## 3. A discordant variant, read mechanistically

**CYP2C19\*2 — rs4244285, chr10:g.94781859G>A.** Rule tripped: **`flag_common_but_plp`** (gnomAD AF above 1% with a drug-response label).

| | |
|---|---|
| Consequence | `synonymous_variant` |
| CADD | **5.69** — bottom band, LR 0.07, *moderate evidence towards benign* |
| REVEL | blank (missense-only) |
| AlphaMissense | blank (missense-only) |
| gnomAD AF | **0.168** |
| ClinVar | drug response, **0 stars** |

This is the variant my section 1 calibration would have called moderately benign. It is also the most important no-function allele in the gene and the reason CYP2C19 genotyping exists in cardiology.

**What would have to be true for both outputs to make sense.** They are answering different questions, and on its own question CADD is right. \*2 is c.681G>A, p.Pro227Pro. The amino acid does not change, so a tool scoring protein-coding constraint correctly reports a low number, and two tools scoring amino-acid substitutions correctly report nothing at all.

The damage is post-transcriptional. c.681G>A creates a **cryptic splice donor site inside exon 5**. The spliceosome prefers it to the authentic downstream donor, which produces a transcript with a 40 bp deletion, a frameshift and a premature stop codon, and therefore a truncated protein with no catalytic activity. The allele is null not because the codon changed but because that stretch of sequence stopped being in the mature mRNA.

**Why I argue for this over the alternatives.** Three other explanations are available and all three lose:

- *ClinVar is wrong and this is a benign polymorphism at 17%.* No. It has CPIC practice-guideline status and drives a clopidogrel boxed warning. The 0-star review status reflects how pharmacogenomic variants get deposited, not weak evidence, and that is a fact about ClinVar's submission conventions rather than about this variant.
- *CADD is wrong.* No, and this is the more useful rejection. CADD is not wrong, it is **out of scope**. Treating this as tool error points at the wrong fix, which would be recalibrating CADD. The right fix is adding a splicing predictor, because no amount of recalibration gives a coding-constraint score the ability to see a splice site.
- *Synonymous variants can change codon usage and translation efficiency.* A real mechanism, but the wrong one here. It predicts a graded, modest effect on protein yield. \*2 behaves as a clean loss of function, which is what a truncating event predicts and what codon-usage effects do not.

Only the splicing explanation predicts both outputs at once: a low protein-level score **and** complete loss of function. It is also the Part 1 toy result appearing in real data — `g.39C>T` was synonymous, scored nothing on consequence or conservation, and produced a donor gain of 1.000.

**What I take from it.** `flag_common_but_plp` didn't fire because something broke. It fired where a frequency-based assumption met a mechanism-based truth, and that collision is where the biology was.

---

## 4. What I could not determine

For this question I switch variants, from \*2 to **L131P (rs371309036, chr10:g.94775450T>C)**. \*2's mechanism is already settled, so there is little I *couldn't* determine about it. L131P is where my pipeline actually committed to something: CADD 24.0, REVEL 0.407, AlphaMissense 0.927, AF 6.6×10⁻⁶, exported as **VUS, leaning pathogenic** on a band LR of 2.14. It also tripped `flag_revel_vs_am`, because the two missense tools differ by 0.52 on it.

**1. Mechanism and direction.** Every score I have is a magnitude with no sign. A proline substitution would be expected to disrupt local backbone geometry, but nothing in my output separates reduced catalytic activity from impaired haem binding from a protein that never folds — or, the case that would reverse the clinical meaning entirely, *increased* activity. \*17 is the standing proof in my own gene that a pharmacogene variant can go the other way, and my calibration labelled \*17 benign. A single scalar from a tool trained on "damaging versus tolerated" cannot tell me which direction the enzyme moved, and for drug dosing the direction is the entire question.

**2. Penetrance, compounded by zygosity and phase.** At AF 6.6×10⁻⁶ this variant has effectively never been seen in a phenotyped carrier, so I cannot say what fraction of carriers show any measurable metabolic difference. CYP2C19 phenotype is also assigned from the **diplotype**, not the variant, so the call is incomplete without knowing whether L131P is heterozygous or homozygous, and if a carrier also has \*17, whether the two are *in cis* (one allele carrying both, the other normal) or *in trans* (one reduced, one increased). Short-read data of the kind behind these annotations is phase-blind. Same variant, opposite predicted phenotype, and nothing in my analysis can choose between them.

*A third worth naming because it is specific to pharmacogenes: **phenoconversion**. A genotypic normal metabolizer taking omeprazole is a phenotypic intermediate metabolizer, because omeprazole inhibits CYP2C19. Genotype is not phenotype when a co-prescribed drug does the same thing the variant does.*

**The experiment that would resolve it.** A **heterologous expression functional assay**. Express recombinant CYP2C19 L131P alongside \*1 wild type and a \*2 null control, then measure activity against a probe substrate — S-mephenytoin 4′-hydroxylation, or clopidogrel bioactivation if clinical relevance is the target. Report Vmax and Km rather than single-point turnover, so that reduced affinity is distinguishable from reduced capacity, and pair it with an immunoblot to separate "folds but is inactive" from "never accumulates". This answers mechanism and direction together, which is exactly what no predictor can do.

**What would change my mind.** Activity at **80% of wild type or above** on both kinetic parameters would falsify my leaning-pathogenic call outright. I would revert L131P to a plain VUS and conclude that AlphaMissense's 0.927 is the false positive in the `flag_revel_vs_am` disagreement, with REVEL's 0.407 the better-calibrated of the two for this gene. Activity **below 20%** would move it from supporting towards functional PS3-level evidence, a strength no amount of extra in-silico scoring could ever reach. Either way the decisive thing is that it is a measurement rather than a prediction, and the distance between my LR of 2.14 and a real number is what this module was about.
