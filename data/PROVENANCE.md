# Where the colon data came from

Everything in this folder was made by `prepare_coad_data.py` in the repo root,
from real downloads. Nothing here is made up or simulated. Running that script
again rebuilds all of it.

## data/coad_transcriptomics.tsv
* **Study** TCGA Colon Adenocarcinoma (TCGA-COAD), TCGA Network, *Nature*
  2012;487:330-337.
* **Files** `TCGA-COAD.star_tpm.tsv.gz` and
  `gencode.v36.annotation.gtf.gene.probemap`, from the UCSC Xena GDC hub at
  `https://gdc-hub.s3.us-east-1.amazonaws.com/download/`.
* **What the source contains** log2(TPM+1) values, STAR quantification, GENCODE
  v36, rows labelled with Ensembl gene IDs.
* **Samples** 471 primary tumours (barcode type `-01`) and 41
  adjacent normals (`-11`). I left out the one recurrent (`-02`) and one
  metastatic (`-06`) sample.
* **What I did** mapped Ensembl IDs to gene symbols, keeping the
  highest-expressed ID where a symbol had several; dropped genes with mean
  log2(TPM+1) below 1 in both groups; took `log2fc` = mean(tumour) -
  mean(normal), which is a real log2 fold change because the values are already
  logged; and used a two-sided Mann-Whitney U test for `pval`, since these are
  different patients and so not paired.
* **Rows** 17611

## data/coad_proteomics.tsv
* **Study** CPTAC-2 Prospective Colon Cancer, Vasaikar et al., *Cell*
  2019;177:1035-1049.e19 (PMID 31031003, doi 10.1016/j.cell.2019.03.030).
* **Files** `...PNNL_Tumor_TMT_UnsharedLogRatio.cct` and
  `...PNNL_Normal_TMT_UnsharedLogRatio.cct`, from the LinkedOmics CPTAC-COAD
  cohort.
* **What the source contains** gene-level TMT log2 ratios against the reference
  channel, using unshared peptides only. Rows are gene symbols.
* **Samples** 96 patients with both a tumour and an adjacent-normal
  channel, so this comparison is properly paired.
* **What I did** took tumour minus normal within each patient, used the mean of
  those differences as `log2fc`, and a Wilcoxon signed-rank test on them for
  `pval`. Genes with fewer than 20 usable pairs were dropped.
* **Check** my log2FC matches the study's own `Tumor_Normal_log2FC.cct` at
  Pearson r = 1.0000 across 6420 shared genes. That's how I know
  this layer is right.
* **Rows** 6709

## data/coad_genomics.tsv
* **Study** the whole-exome data from the same CPTAC-2 colon study, so the same
  patients as the protein layer. WUSM GATK pipeline, gene-level calls.
* **File** `...WUSM__Mutation__GAIIx__03_01_2017__BCM__Gene__GATK_Pipeline.cbt`
  from the LinkedOmics CPTAC-COAD cohort.
* **What the source contains** a 1 if that patient's tumour has a somatic
  mutation in that gene, 0 if not. 106 patients.
* **Gene lengths** longest CDS per gene symbol, from UCSC hg38 `refFlat.txt.gz`.
* **What I did** `mut_freq` is the fraction of patients mutated. Then I worked
  out a background rate of mu = 1.602e-05 mutated-gene events per coding base
  per patient across all genes with a known CDS length, converted it to a
  per-gene chance of being hit as p = 1 - exp(-mu x CDS length), and used a
  one-sided binomial test on the observed number of mutated patients.
  `neglog10p` is -log10 of that p-value.
* **Why bother with the length correction** without it the top of the list is
  TTN (56.6% of patients), MUC16 (36.8%) and OBSCN (27.4%), all of which
  outrank or match KRAS (33.0%) just by being long. This is a length
  correction in the same spirit as MutSigCV, not a proper driver-gene method,
  and I don't present it as one.
* **Rows** 14706

## Quick check that the mutation file is what it says it is
Straight out of the CPTAC file, the gene-level frequencies are APC 77.4%,
TP53 52.8%, KRAS 33.0%, BRAF 17.0%. That's the standard colon adenocarcinoma
picture, which is reassuring.

## data/coad_cptac_matched_rna_protein_corr.tsv (context, not used for scoring)
Per-gene Spearman correlation between CPTAC tumour RNA
(`...UNC__RNAseq...BCM_RSEM_UpperQuartile_log2.cct`) and CPTAC tumour protein,
across the 96 tumours measured on both. 7069 genes
with at least 20 usable samples, median Spearman r = 0.336. This is
only here so I can say what the cross-cohort design costs, using numbers from
inside one cohort.

## data/coad_kras_variant_spectrum.tsv (context, not used for scoring)
From the site-level calls in
`...WUSM__Mutation__GAIIx__03_01_2017__BCM__Site__GATK_Pipeline.cbt`.
35 of 106 patients have a KRAS coding mutation, across
11 different variants.

## data/coad_independent_proteomics_fudan.tsv (context, not used for scoring)
Supplementary Table S6B of "Comprehensive Proteogenomic Profiling Reveals the
Molecular Characteristics of Colorectal Cancer at Distinct Stages of
Progression", *Cancer Research* 2024;84:2888-2910. A separate Chinese (Fudan)
cohort with paired tumour and adjacent-normal proteomics, measured label-free
instead of by TMT. 28 complete pairs, 4171 proteins
with at least 5 usable pairs. Missing values are written as 1e-05 in that
table and I treat them as not detected. I use it to check the CPTAC protein
layer against a completely separate group of patients.

## About the two different cohorts
The RNA layer is TCGA-COAD; the protein and mutation layers are CPTAC-2 colon.
**No patient is in both.** So everything is joined at the gene level and I
can't say anything about individual patients. CPTAC does have RNA-seq on its
own 106 patients, but only for the tumours, with no adjacent-normal RNA, so it
can't give a tumour vs normal RNA comparison. That's the reason the RNA layer
had to come from TCGA.
