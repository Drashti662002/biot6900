# PRS portability for colorectal cancer across 1000 Genomes ancestries

**BIOT 6900 · Module 6 Assignment · Drashti Bharat Bhanushali · 29 October 2026**

*I carried my disease forward from Modules 2–4. In those modules it was colon adenocarcinoma, and PGS Catalog scores are published for colorectal cancer as a whole, so that is the trait here. Module 2 looked at somatic changes inside tumours (KRAS, TP53, APC). A PRS is a different question: it uses common inherited variants to estimate who is more likely to develop the cancer at all.*

*The score is **PGS000148** (CRC63, Jeon et al., Gastroenterology 2018). It has 63 SNPs, with per-allele log-odds weights estimated in **European-ancestry samples only**. Canvas `data.zip` holds only the 22 Alzheimer's SNPs, so I rebuilt the four input files myself (`prepare_crc_data.py`) from the same sources the professor used: the PGS Catalog harmonized file and 1000 Genomes phase 3 (GRCh37). I used the same 2,336 people, ancestry labels and ages as the class cohort. **The phenotype is simulated, and so is age.** I used the professor's recipe from the Canvas README with every number unchanged (seed 6900, ±250 kb, r² > 0.6 in EUR, noise SD 2.09, top 30% are cases), and I did not retune anything. The PRS and its distribution across ancestries come from real genotypes and real weights. Every AUC and calibration number depends on the simulated label, so those numbers illustrate the effect but do not measure real colorectal cancer risk.*

---

## Proof of data (Section 2 output, pasted verbatim)

```
================================================================
  RUNNING ON REAL DATA
================================================================
source      : REAL
people      : 2336
score size  : 62 variants
by ancestry : {'AFR': np.int64(500), 'EAS': np.int64(500), 'EUR': np.int64(500), 'SAS': np.int64(489), 'AMR': np.int64(347)}
  genotypes.csv: sha256 5617fc3c
  variants.csv: sha256 0b8c8747
  score.csv: sha256 309987fb
  meta.csv: sha256 8dbbb462
```

`USE_DEMO = False`. The checksums match `data/README.md`, which `prepare_crc_data.py` writes. They differ from the Canvas README's checksums because these are colorectal files, not the Alzheimer's ones.

## Results at a glance

**Harmonization (3.1): used 61, flipped 28, dropped 1.** Of the 63 published SNPs:
- **rs202110856 has no record in 1000 Genomes phase 3** at its GRCh37 position, so it can't be scored.
- **rs7229639 (SMAD7) is dropped by the keep/flip/drop rule.** The published row has a blank `effect_allele`, so the score never says which allele the 0.023 weight belongs to. It matches neither REF nor ALT, so the rule drops it. I kept that locus in the *simulated* disease, so the score loses a real (small) signal, as it would in practice.
- **Three SNPs are strand-ambiguous** (A/T or C/G: rs812481, rs35360328, rs16941835). The keep/flip/drop rule can't detect a strand flip at these. I kept them because both the PGS Catalog file and the VCF are on the + strand.
- **Alignment check:** EUR cases score higher than EUR controls (AUC 0.685), so the effect alleles are lined up correctly.

| | EUR | AFR | EAS | SAS | AMR | Overall |
|---|---|---|---|---|---|---|
| n | 500 | 500 | 500 | 489 | 347 | 2,336 |
| Mean PRS (SD units) | +0.16 | −0.09 | −0.17 | **+0.30** | −0.28 | 0 |
| % above pooled top-10% cutoff (PRS > 1.282) | 12.2 | 8.4 | 7.2 | **14.9** | 6.3 | 10.0 |
| % cases (simulated) | 30.0 | 26.6 | 28.2 | 38.0 | 26.2 | 30.0 |
| AUC (3.4) | **0.685** | **0.567** | 0.689 | 0.606 | 0.670 | 0.647 |
| AUC 95% bootstrap CI | 0.633–0.735 | 0.511–0.625 | 0.638–0.737 | 0.555–0.656 | 0.602–0.732 | |
| Mean tag–causal r² | 0.83 | **0.45** | 0.76 | 0.77 | 0.74 | |
| Observed / predicted (EUR-fitted, 3.5) | 1.00 | 0.99 | 1.10 | **1.20** | 1.08 | |

Incremental value (3.6, 5-fold CV, pooled): **age only AUC 0.579, age + PRS 0.666** (+0.087).
For comparison, the published EUR validation AUROC for PGS000148 is 0.59 (PGS Catalog PPM000450, with age, family history and endoscopy as covariates).
Figures: `figures/prs_by_ancestry.png`, `figures/calibration_by_ancestry.png`. Every extra number below comes from section 3.8 of the notebook.

---

## 1. Why the European-derived score loses accuracy in other ancestries

By construction, the simulation gives each causal SNP **the same effect in every ancestry**, and one liability rule applies to everyone. Biology is held constant, so any difference between groups comes from LD and allele frequency.

**LD: the score SNPs are tags, and tags only work where LD holds.** The 63 SNPs are GWAS hits chosen because they correlate with a causal variant in Europeans. A tag carries signal in proportion to its r² with that causal variant *in the population it's used in*.

| | EUR | AFR | EAS | SAS | AMR |
|---|---|---|---|---|---|
| Mean tag–causal r² | 0.83 | **0.45** | 0.76 | 0.77 | 0.74 |
| Tags with r² < 0.2 | 0 | **25 of 62** | 7 | 2 | 2 |

African populations have the largest long-term effective population size, so recombination has broken haplotypes into shorter LD blocks. A tag that sits on the causal variant's haplotype in Europe often doesn't in Africa. For 25 of the 62 tags, AFR r² falls below 0.2. AFR is also the only group whose AUC interval lies entirely below EUR's (0.511–0.625 vs 0.633–0.735).

**Allele frequency decides whether an LD loss costs anything.** Several of the largest-weight tags lose LD in more than one ancestry, but the damage depends on how common the tag is there:

| Tag (weight) | r² EUR / AFR / EAS | Effect-allele freq EUR / AFR / EAS |
|---|---|---|
| rs16892766, EIF3H 8q23 (0.200) | 0.69 / **0.06** / 0.01 | 0.09 / **0.14** / 0.001 |
| rs73208120, 12q24 (0.188) | 0.64 / 0.16 / undefined | 0.085 / 0.005 / 0.000 |
| rs1535, FADS1/2 (0.184) | 0.97 / **0.12** / 0.99 | 0.65 / 0.89 / 0.44 |

In EAS, rs16892766 and rs73208120 are almost monomorphic. A SNP that doesn't vary adds nothing to the score, so poor tagging costs nothing either. In AFR, rs16892766 is *more* common than in Europe (14% vs 9%) but tracks its causal SNP at r² 0.06. That puts one of the largest weights in the score on a variable that is essentially noise. The same r² loss is harmless in one population and harmful in another, depending on allele frequency.

**Frequency also explains why EAS matches EUR here.** The largest weight, rs11064437 (0.218), had no LD neighbour, so the simulation made it its own causal SNP, with r² = 1 everywhere by design. Its effect allele is almost fixed in EUR (99%), so it contributes almost nothing there. In EAS it sits at 72% and varies, so it adds real, perfectly tagged signal. This is partly an artifact of the simulation: in reality the causal variant would rarely be the GWAS lead SNP itself. I would not read "EAS = EUR" as the score being portable to East Asians. EUR and EAS AUCs are statistically indistinguishable here, and that is the most the data supports.

**SAS is the case I can't fully explain.** Its mean r² (0.77) is as good as EAS's, but its AUC is lower (0.606, CI 0.555–0.656). The intervals overlap, so part of this is noise. Mean r² also weights every tag equally, while the AUC depends on which tags carry variance in that group. That points to a limit of the analysis: a single averaged r² is not enough to predict portability.

**Allele frequency also moves the whole distribution.** The PRS adds up effect alleles, so a group that carries more of them scores higher whether or not its risk is higher. The net shifts here are modest (EUR +0.16, SAS +0.30, AFR −0.09, EAS −0.17, AMR −0.28), but they hide large opposing per-SNP effects. In AFR, rs11064437 alone moves the mean by −0.41 SD and rs4246215 (FADS) by −0.35, while rs1535 (+0.27), rs174550 (+0.27), rs4779584 near GREM1 (+0.26) and rs3184504 (+0.24) push the other way. A group's average PRS is therefore partly an accident of how these cancel. A score built from joint-model weights with opposite signs inside one locus, as at FADS1/2, is especially sensitive to this.

Two more structural reasons: the SNPs were selected in a European GWAS, so loci that matter mainly in other ancestries were never found. And every SNP the score misses (rs202110856, rs7229639) reduces accuracy everywhere.

## 2. One threshold, many groups

**What the pooled top-10% cutoff (PRS > 1.282) does:**

| | EUR | AFR | EAS | SAS | AMR |
|---|---|---|---|---|---|
| % of group flagged | 12.2 | 8.4 | 7.2 | 14.9 | 6.3 |
| Share of all 233 flags | 26.1% | 17.9% | 15.4% | 31.2% | 9.4% |
| Case rate if flagged | 57.4% | **26.2%** | 55.6% | 57.5% | 54.5% |
| Case rate if not flagged | 26.2% | **26.6%** | 26.1% | 34.6% | 24.3% |
| Share of the group's cases caught | 23.3% | **8.3%** | 14.2% | 22.6% | 13.2% |
| Group's own top-10% PRS | 1.38 | 1.17 | 1.10 | 1.51 | 1.05 |

**It isn't fair, but not simply because flag rates differ.** Different flag rates can be acceptable if the groups really differ in risk. In this simulation SAS really has more cases (38%), so flagging more SAS people is partly correct. The better fairness test is whether a flag means the same thing in every group. For EUR, EAS, SAS and AMR it roughly does: a flagged person's case rate is 55–58% in all four. **For AFR the flag means nothing.** A flagged African-ancestry person has a 26.2% case rate against 26.6% if not flagged. The cutoff still fires for 42 AFR people, and every one of those flags is uninformative.

The second unfairness is sensitivity. The cutoff catches 23% of EUR cases, 14% of EAS, 13% of AMR and **8% of AFR**. Even where a flag is meaningful, it reaches fewer of the people who will actually get the disease.

**If a clinic applied it unchanged**, the most likely real use for a CRC score would be risk-adapted screening: earlier or more intensive colonoscopy for people above the cutoff. That was the stated purpose of Jeon et al.'s paper.

- **AFR patients would be triaged at random.** 42 would get earlier colonoscopy for no reason, and 92% of AFR cases would get no earlier screening. That matters more for this disease because, in the US, Black Americans have among the highest colorectal cancer incidence and mortality of any group. An unchanged EUR score could put the highest-burden group at the back of the queue.
- **EAS and AMR** flags are informative but rare, so most of their cases would be missed (86–87%).
- **The more dangerous version is the reverse use**: treating a *low* score as a reason to delay screening. That turns a missed flag into a missed cancer.

**Group-specific cutoffs** (last row) would equalize flag rates, but they can't add accuracy. In AFR they would still flag at random. Super-population labels are also crude: AMR and AFR include admixed populations (MXL, PUR, ASW, ACB), and many patients don't fit any one category. Standardizing the PRS against continuous genetic ancestry (principal components), as the PGS Catalog calculator and eMERGE do, handles the location shift. **It does not fix the AFR accuracy problem. Equal flag rates are not the same thing as equal accuracy.**

## 3. Discrimination vs calibration

**Discrimination (3.4)** asks whether the score ranks a random case above a random control *within* a group. Shifting the whole distribution up or down doesn't change it. Here, the score ranks moderately well in EUR, EAS and AMR (0.67–0.69), less well in SAS (0.606), and only a little better than chance in AFR (0.567). AUC picks up the LD problem but ignores the frequency problem.

**Calibration (3.5)** asks whether a predicted 30% means 30%. I fitted the logistic model in EUR only (intercept −1.054, slope 0.720 per SD) and applied it to everyone.
- **EUR:** looks perfect (0.300 vs 0.300), but by construction, since the model was fitted on those same people. EUR is the training set, not a validation.
- **AFR:** right on average (0.268 predicted vs 0.266 observed), but **flat**. Predicted risk rises from 0.11 to 0.47 across quintiles while observed risk only moves between 0.16 and 0.34.
- **SAS:** under-predicted overall (0.318 vs 0.380, O/E 1.20), and worst in the bottom two quintiles: predicted 0.15 and 0.23, observed 0.32 and 0.32. SAS people with low scores have about twice the risk the model gives them.
- **EAS and AMR:** close, slightly under-predicted (O/E 1.10 and 1.08), and their curves follow the diagonal.

**Where they disagreed:**

1. **AFR passes calibration-in-the-large and fails everything else.** An O/E of 0.99 is the best of any group except the training set. Judged on mean calibration alone, AFR would be the most "honest" group. The AUC (0.567) and the flat calibration slope show the model has almost no ranking ability there. The average is right only because over- and under-predictions cancel. At the top end, flagged AFR people are predicted at **55.0%** risk and observed at **26.2%**. For flagged EUR people the figures are 55.7% and 57.4%. The model gives both groups the same number, and it is wrong by a factor of two in one of them.
2. **EAS ranks as well as EUR but is mis-calibrated.** Same AUC (0.689 vs 0.685), yet 10% under-predicted overall. Good ranking does not guarantee correct absolute risks. EAS mostly needs a new intercept.
3. **SAS has the worst calibration but only a moderate AUC drop.** Its baseline risk differs from EUR's, and a model built on the PRS alone can't learn that. SAS needs both a new intercept and a new slope. AFR needs a better score.

Caveat: each calibration point is a quintile of about 70–100 people, so each observed rate has a standard error of about ±5 percentage points. The AFR flatness, the AFR top-end gap and the SAS low-end gap are well outside that. The small wiggles in the EUR and AMR curves are not.

## 4. Clinical utility and limits

**What the score adds to age.** Pooled, age alone gives a cross-validated AUC of 0.579 and age plus PRS gives 0.666, an increase of 0.087. Within each ancestry the gain is uneven:

| | EUR | AFR | EAS | SAS | AMR |
|---|---|---|---|---|---|
| Age only | 0.604 | 0.563 | 0.452 | 0.649 | 0.553 |
| Age + PRS | 0.708 | 0.581 | 0.682 | 0.680 | 0.679 |
| Gain | +0.104 | **+0.018** | +0.229 | +0.031 | +0.126 |

In AFR the score adds almost nothing on top of age. EAS's age-only AUC of 0.452 is below chance. Age has the same simulated effect in every group, so that number is sampling noise in a 500-person group, and EAS's large gain should be read with that in mind. Both the age effect (0.35 per 8 years) and the phenotype are simulated, so none of this measures real incremental value. Real CRC risk rises much more steeply with age than this.

**Where it could help:**
- **Risk-adapted screening in European-ancestry adults**, the use the score was built for: moving the start age of colonoscopy or FIT earlier for people in the top of the distribution. CRC is a good fit because there is an effective, cheap intervention (screening and polyp removal), so a modest AUC can still be useful when it changes *when* someone is screened rather than *whether* they're treated.
- **Choosing between screening methods**, for example colonoscopy for high-PRS people and FIT for average-risk people, in a validated population.
- **Research:** enriching prevention trials (aspirin, diet) for higher-risk European-ancestry participants.

**Where it should not be used:**
- **To delay or skip screening for a low score**, in any ancestry. Most CRC cases occur in people with average scores, and guideline screening from age 45 already covers everyone.
- **In African-ancestry patients**, where it barely discriminates, and in admixed patients, where it hasn't been validated.
- **As a substitute for hereditary cancer evaluation.** Lynch syndrome (MLH1, MSH2, MSH6, PMS2) and FAP (APC) are rare high-penetrance variants that a 63-SNP common-variant score can't detect. Family history or early-onset cancer calls for germline panel testing, not a PRS.
- **To say anything about a tumour.** Module 2's KRAS, TP53 and APC findings are somatic, and this score is germline.
- **For insurance or employment decisions**, or direct-to-consumer reporting without counselling.

**What I would need before deploying it:**
1. **A real phenotype.** Everything in 3.4–3.6 rests on a simulated label and a simulated age. Validation needs registry- or pathology-confirmed incident colorectal cancer, with age at diagnosis, sex, family history, screening/colonoscopy history (polyp removal prevents cases, which biases case ascertainment) and subsite (colon vs rectum, which matters for a colon adenocarcinoma question), in genotyped people with enough cases per ancestry. Candidate datasets:
   - **UK Biobank** and **FinnGen** (registry-linked; mostly European).
   - The **GECCO / CORECT / CCFR** consortia (the CCFR includes family-history and multi-ethnic sites).
   - The **Multiethnic Cohort** (African American, Japanese American, Latino and Native Hawaiian participants).
   - The **Southern Community Cohort Study** (mostly African American).
   - **All of Us** and the **Million Veteran Program** (EHR outcomes across ancestries).
   - **Biobank Japan**, **China Kadoorie Biobank** and the **Asia Colorectal Cancer Consortium**.
   - African-ancestry CRC GWAS collections. These are still small, and that is the main gap.
2. **Multi-ancestry weights.** The PGS Catalog already has candidates: Fernandez-Rozadilla 2022 (PGS003850, EUR + EAS) and Thomas 2023 (PGS003852, LDpred2 on EUR + EAS). Neither has enough African-ancestry training data, so new GWAS is needed. Methods like PRS-CSx can combine ancestries.
3. **Calibration within each ancestry** (intercept and slope) and continuous-ancestry standardization, judged by calibration slope and by the flagged-group case rate, not only overall O/E.
4. **Absolute 10-year risk at a stated age**, accounting for death from other causes, combined with family history and lifestyle as Jeon et al. did. "Top 10%" is not something a patient can act on.
5. **Prospective evidence** that PRS-guided screening finds more early cancers or advanced adenomas per colonoscopy than age-based screening, without widening disparities.
6. **Governance**: consent, counselling, data protection, and explicit rules against using a low score to deny or delay screening.

---

## Note on the simulation (artifact vs biology)

Some patterns here come from how the simulation was built, not from biology:
- **EUR LD advantage built in:** causal SNPs were chosen from SNPs in LD with each tag *in Europeans*, so the score tags its causal SNPs best in EUR by construction. The AFR drop uses real 1000 Genomes LD, but its size depends on this design.
- **Own-causal tags transfer perfectly:** five tags (including the largest weight, rs11064437) had no qualifying neighbour and became their own causal SNPs. They transfer perfectly to every ancestry, which flatters portability (especially for EAS).
- **Case rates reflect simulated frequencies:** SAS's higher case rate (38%) comes from its causal-allele frequencies in this simulation. It is not a claim about real colorectal cancer incidence in South Asians, which is lower than in Europeans.

None of this is evidence about how colorectal cancer biology differs between ancestries. In this dataset, by design, it doesn't.

**Ethics.** Almost all the GWAS data behind this score comes from European-ancestry participants. Used unchanged in a mixed clinic, it would give the most accurate screening decisions to the group with the most research behind it, and essentially random ones to one of the groups with the highest colorectal cancer mortality in the US. The fix is mainly better data: diverse GWAS cohorts, validated by ancestry, recruited with community involvement. Adjusting thresholds afterwards can't replace that.
