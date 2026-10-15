# Target triage with protein language models in colon adenocarcinoma

**BIOT 6900 · Module 3 Assignment · Drashti Bharat Bhanushali · 15 October 2026**

*I carried my Module 2 disease forward, so these are colon adenocarcinoma targets rather than Alzheimer's ones. The method is exactly the one the lab asks for.*

I took the top 15 genes from my Module 2 multi-omics ranking (`targets_coad.csv`, 4,851 genes scored under equal weights), looked each one up in the data pack proteome, and scored them with a tractability prior built from ESM C embeddings and Open Targets clinical precedence labels. All 15 genes mapped to a reviewed protein, so I have no unmapped genes to report. Only one of them, CA2, was in the 1,201-protein training panel, and I removed it before fitting anything, so nothing below was predicted by a model that had seen its own label.

---

## 1. Re-ranking: did the prior agree with Module 2?

**No, and the numbers say so clearly.** The Spearman correlation between my Module 2 percentile and the tractability percentile is **−0.089 (p = 0.752, n = 15)**. That is indistinguishable from no relationship. The two scores are measuring genuinely different things, which is the whole reason the re-ranking step exists. If they had agreed, the prior would have added nothing.

The biggest riser was **TSPAN7**, which went from 9th to 4th. The biggest fallers were **DES** (2nd to 6th) and **PRPH** (8th to 12th), both down four places.

DES and PRPH are worth looking at together because they fell for the same reason. Both are in the **intermediate filament family** (DES is P17661, PRPH is P41219), and both scored **0.001 on both priors** — the lowest values anywhere in my table. Intermediate filaments are structural fibres. They polymerise into long coiled-coil ropes, they have no catalytic site, and nothing about them looks like the proteins that people have historically managed to drug. The model has never been told any of that. It worked it out from sequence resemblance alone, and it came to the same conclusion a medicinal chemist would. That is the single most convincing thing the prior did on my list.

CA2 stayed at rank 1. I checked whether that was an artefact of my weight by recomputing the ranking across the full range of `w`, and CA2 is first at every value from 0 to 1. It is first on Module 2 evidence and first on tractability.

### Justifying w = 0.6

I used **w = 0.6**, so Module 2 evidence carries slightly more weight than the prior.

My reasoning is that the two inputs are not the same kind of evidence. The Module 2 score comes from three independent measurement layers in real patients — 471 tumour versus 41 normal RNA samples, 96 matched protein pairs, and 106 exomes. The tractability prior is one model's inference about sequence resemblance to historically drugged proteins, and I can state exactly how well it does: family-split AUROC 0.751. Measurement should outweigh inference, but not so heavily that the prior cannot move anything.

I also checked how much the choice actually matters, because I did not want to assert a number I had not tested:

| w | Top three | TSPAN7 | DES | Max shift vs w = 0.6 |
|---|---|---|---|---|
| 0.0 | CA2, TSPAN7, BCHE | 2 | 14 | 8 places |
| 0.4 | CA2, BCHE, TSPAN7 | 3 | 9 | 3 places |
| **0.6** | **CA2, BCHE, CFL2** | **4** | **6** | **—** |
| 0.8 | CA2, DES, BCHE | 7 | 2 | 4 places |
| 1.0 | CA2, DES, PTN | 9 | 2 | 5 places |

A gene can move eight places across the range, so `w` is a real decision and not a cosmetic one. I should be honest that any value between about 0.4 and 0.7 would be defensible on the argument I just made, and I picked the middle of that range. What does not change is CA2 at the top, and what changes most is DES, which is exactly the gene the two scores disagree about most.

---

## 2. Modality: CR2 is the one clear split

**CR2** is the only target where the antibody prior beats the small-molecule prior by a meaningful margin: **p_ab = 0.663 against p_sm = 0.142**, a gap of +0.521. The next largest gap in the same direction is COLEC12 at +0.121, so CR2 really is on its own.

The biology lines up with the numbers. CR2 (P20023) is a **1,033-residue cell-surface receptor** in the receptors of complement activation family. It is a long extracellular protein sitting on the outside of the cell, built from repeated domains, with no deep enclosed cavity for a small molecule to occupy. But the whole of it is exposed to the outside of the cell, which is precisely what an antibody needs. Its three nearest drugged neighbours in the panel are **CSF2RA (0.94), IFNGR2 (0.94) and IL1RL2 (0.94)** — all cytokine receptors, all surface proteins, all the kind of thing antibodies are raised against.

CA2 is the clean opposite case, and the contrast makes the point better than CR2 alone: **p_sm = 0.770, p_ab = 0.008**, a gap of −0.761. CA2 (P00918) is a 260-residue soluble enzyme with a deep zinc-containing active site, and its nearest neighbours are **CA1 (0.99) and CA7 (0.96)**, other carbonic anhydrases. Carbonic anhydrase inhibitors like acetazolamide are approved drugs, so this is a family with real small-molecule precedent, and the model found it.

**One caution I want to record.** In the training panel, only **5.6% of proteins are antibody-positive** against **25.2% small-molecule-positive**. The antibody prior is fitted on roughly a fifth as many positive examples, so it is the less reliable of the two. I would treat CR2's 0.663 as a flag worth following up rather than a confident probability.

---

## 3. Validation: the honest number is the family split

| Split | AUROC | Average precision |
|---|---|---|
| Random (StratifiedKFold) | 0.833 | 0.707 |
| **Family (GroupKFold)** | **0.751** | **0.578** |

Chance level for average precision is the positive fraction, **0.333**.

The random split is dishonest here for a specific reason. Protein families are full of close relatives — the panel contains CA1, CA2 and CA7, and several EPH receptors. Shuffle proteins at random and a carbonic anhydrase in the test fold will have two of its siblings sitting in the training fold. The model does not need to understand tractability to score well; it only needs to recognise a family it has already been shown. The family split removes that crutch by holding out whole family groups, and performance drops by **0.082 AUROC and 0.129 average precision**.

The family-split number is the one I would give a discovery team, because the question they are actually asking is "will this work on a protein class we have not characterised yet?", and that is the question the family split answers.

I would add one caveat that makes even 0.751 look slightly optimistic. The panel has **821 family groups across 1,201 proteins, and 276 of those proteins are singletons** with no family annotation. A singleton is its own group, so holding it out is no different from a random hold-out. For roughly 23% of the panel the family split provides no protection at all, which means the true out-of-family performance is probably a little below 0.751 rather than above it.

I report average precision alongside AUROC because the panel is enriched: a third of it is positive, far above the natural rate across the proteome. Accuracy would look impressive and mean nothing. Average precision has a chance level you can state (0.333), so 0.578 is interpretable — better than chance, not dramatically so.

---

## 4. Anchor variants

All four scores come from ESM-2 650M as `log P(mutant) − log P(wild-type)`, in UniProt precursor numbering. The percentile is relative to every possible substitution in that same protein.

| Variant | Note | Score | Percentile |
|---|---|---|---|
| APOE R176C | rs7412, defines ε2 (R158C mature) | −8.56 | 11.8 |
| APOE C130R | rs429358, defines ε4 (C112R mature) | **+7.21** | **100.0** |
| TREM2 R47H | rs75932628, AD risk variant | −5.16 | 47.7 |
| TP53 Y220C | hotspot, druggable crevice | −3.80 | 33.4 |

**The APOE pair is the most instructive result in this whole assignment.** C130R defines the ε4 allele, which is the best-established common genetic risk factor for late-onset Alzheimer's disease. The model gives it a score of **+7.21, at the 100th percentile** — meaning that of every possible substitution in APOE, this is the one the model finds *most* likely. It does not merely fail to flag ε4; it actively prefers it. The reason is that ε4 is the ancestral allele shared with other primates, so it is the residue the evolutionary record supports. Meanwhile ε2 (R176C), which is generally protective, scores −8.56 at the 11.8th percentile, as though it were damaging.

That is a complete inversion of clinical risk, and it is not a bug. It is the model doing exactly what it was trained to do.

**Why a strongly negative score is not the same as pathogenic.** The score answers one question: is this residue the kind of thing evolution kept at this position? Evolution selects on reproductive fitness across deep time, not on human disease. A variant that causes late-onset disease, after reproduction, is nearly invisible to selection. Gain-of-function variants often sit at poorly conserved positions, because they create new activity rather than destroying old activity, so the model is not surprised by them at all. And conservation can reflect constraints that have nothing to do with the phenotype in question. TREM2 R47H makes the point from the other side: it is a confirmed AD risk variant and it scores at the **47.7th percentile**, the middle of the distribution. If I had set a threshold anywhere sensible, I would have missed it.

**What these models cannot see.** They score missense changes only, so synonymous, splice-site, promoter, UTR and non-coding variants return nothing at all — not a benign call, nothing. They have no cell type, no expression level, no patient, no zygosity and no structure. And they cannot see pockets that only exist in the mutant: TP53 Y220C scores −3.80 at the 33.4th percentile, the least extreme of all seven p53 hotspots in the data pack, even though the crevice its mutation opens is a well-known drug-discovery opportunity. The model has no way to know that losing a bulky tyrosine leaves a hole worth targeting.

---

## 5. Limitations, and the two targets going to Week 4

**What sequence alone cannot tell me.** It cannot tell me whether a protein has a pocket, because a pocket is a three-dimensional property and the prior only ever saw a string of letters. It cannot tell me whether the protein is expressed in the tissue I care about. It cannot distinguish inhibition from activation. And its labels encode history rather than biology — a protein scores well largely because it resembles proteins that people already succeeded in drugging, which means the prior is systematically conservative about anything genuinely new. In the Part 1 toy model this showed up as the Channel-like family, which is druggable but scored 0.40 because its pocket appeared nowhere in training. In Part 2 it showed up on p53 itself, which scored 0.11 despite having real clinical precedence.

**I am sending CA2 and CR2 to Week 4.**

**CA2** is my top-ranked target on both lines of evidence, has the highest small-molecule prior on my list (0.770), and sits 0.99 cosine from CA1, a protein with clinical precedence. It is also the one target where I can actually check the structural question properly: CA2 is UniProt **P00918**, which is the same protein as the Module 4 demo, where acetazolamide from crystal 3HS4 fits into the AlphaFold-predicted pocket with no clash. So I am not going into Week 4 guessing whether there is a pocket. I can test whether the prior's confidence is justified against a structure I can inspect.

**CR2** goes for the opposite reason: it is the only modality split I found, and a structure-based assessment would ask a different question of it. For CR2 the question is not "is there a pocket" but "is there an accessible, well-ordered surface epitope", which is a question about domain confidence and inter-domain arrangement across 1,033 residues. Given that antibody labels are the sparse ones in this panel, CR2 is the prediction most in need of independent checking.

**The honest caveat on both choices.** My Module 2 report already noted that 43 of my top 50 genes are *down* in tumour on both layers, and the top 15 I carried forward here are all in that category. You cannot inhibit a protein the tumour has already lost. CA2 is down 5.20 log2 at RNA and 2.16 at protein; CR2 is down 1.97 and 1.42. Both would be better framed as diagnostic or stromal-context markers than as inhibition targets.

I ran a supplementary check on this rather than leaving it as an assertion. I re-scored the seven genes from my Module 2 report that are **up** in tumour on both layers (PLAU, TMEM97, SOX9, LRRC15, MMP11, DPEP1, THBS2), refitting both priors with all 22 scored genes held out so the comparison stays leakage-free:

| Direction | n | Mean best-modality prior | Max |
|---|---|---|---|
| Down in tumour | 15 | 0.291 | 0.788 (CA2) |
| Up in tumour | 7 | 0.255 | 0.513 (MMP11) |

The two groups score about the same, which tells me something useful: **the tractability prior is completely blind to direction of change.** It does not know, and cannot know, whether the tumour has more or less of a protein. That is a limitation of the prior, not of my Module 2 ranking, and it means the prior cannot rescue a target list that has a direction problem — it can only tell you which of the genes you handed it look druggable.

If I were choosing targets for a real programme rather than answering this assignment, I would send MMP11 (p_sm 0.513) and TMEM97 (p_sm 0.403) forward instead, because they are at least present in the tumour. I am sending CA2 and CR2 because they are the top of the ranked list the assignment asked me to build, and because CA2 gives me a structure I can genuinely verify against.

---

## Unmapped genes

**None.** All 15 of my top targets had a reviewed protein entry in `proteome.tsv` and mapped cleanly to an embedding row.

Had any failed to map, the likely causes would be non-coding RNAs, pseudogenes, or genes with no reviewed UniProt entry in release 2026_03, and none of those can be scored by a protein language model at all, because there is no protein sequence to embed. For those I would need a different kind of evidence entirely — expression-level and regulatory data rather than protein tractability — and I would report them explicitly rather than let them disappear from the output, since a gene that silently drops out of a target list looks identical to a gene that was never a candidate.
