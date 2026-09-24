# BIOT 6900 coursework — Drashti Bharat Bhanushali

## Module 1
I used Python and Biopython to check my environment, translate a DNA
sequence, and retrieve information from PubMed, UniProt, PDB, and
the GWAS Catalog.

Part D was skipped as instructed by the professor.

## Module 2 — Multi-Omics Target Identification & Validation

**I did Part 3 on colon adenocarcinoma instead of Alzheimer's**, which my
professor approved beforehand. The method is unchanged: my layers come from
different cohorts with no shared patients, so I join them at the gene level,
which is what Part 3 is about.

I did not edit any of the professor's code. `load_ad`, `rank_percentile`,
`multi_evidence_score`, `EQUAL_WEIGHTS` and all of Parts 1 and 2 are exactly as
he wrote them. I only filled in the `# TODO` cells and added my own cells
around them.

### What's here

| File | What it is |
|---|---|
| `BIOT6900_Module2_Starter.ipynb` | The notebook. Runs top to bottom, checked with Restart & Run All. |
| `targets_coad.csv` | The Part 3 output: all 4,851 genes ranked. This is the file the guide calls `targets_ad.csv` — I named it for my disease. |
| `Module2_Part3_Report.md` | The 3.6 write-up. |
| `prepare_coad_data.py` | Downloads the TCGA and CPTAC files and builds the three summary tables. |
| `data/PROVENANCE.md` | Where every number came from: URLs, file names, sample counts, formulas. |

### Where the data comes from (all open access, no data use agreement)

* **RNA** — TCGA-COAD, 471 primary tumours vs 41 adjacent normals. GDC STAR
  log2(TPM+1), downloaded from the UCSC Xena GDC hub.
* **Protein** — CPTAC-2 Prospective Colon Cancer (Vasaikar et al., *Cell*
  2019;177:1035–1049.e19, PMID 31031003), 96 matched tumour/normal pairs.
  PNNL TMT gene-level log ratios, downloaded from LinkedOmics.
* **Mutations** — the exome data from that same CPTAC study, so the same 106
  patients as the protein layer. WUSM GATK gene-level calls, also from
  LinkedOmics.

None of these come with a differential table already made, so I calculated each
comparison myself (Mann-Whitney for RNA, paired Wilcoxon for protein, and a
length-corrected binomial test for the mutations).

I checked the protein layer twice. My log2FC matches CPTAC's own published
`Tumor_Normal_log2FC` file at Pearson r = 1.0000 over 6,420 genes, which shows
I did the arithmetic right. Then I compared against a completely separate
Chinese cohort (supplementary table S6B of *Cancer Research* 2024;84:2888–2910,
28 paired samples, label-free instead of TMT): 65.5% of proteins move the same
direction in both, Pearson r = 0.451. That's moderate, not great, and the
report says so.

### Running it

```bash
conda activate biot6900
python prepare_coad_data.py     # one time, downloads about 175 MB into data/raw/
jupyter notebook                # then Kernel → Restart & Run All
```

`data/raw/` is in `.gitignore` because the TCGA expression matrix alone is
135 MB. Re-running the script downloads it again.

The last step of the script reads a journal supplementary spreadsheet that
can't be downloaded automatically. If you want it, save
`can-23-1878_supplementary_tables_1-5_suppst1-5.xlsx` into `data/raw/` and
install `openpyxl`. If it isn't there the script skips that step and everything
else still works.

### What didn't work

* **I couldn't get all three layers from one cohort.** I wanted everything from
  CPTAC-2 colon since those are the same patients. CPTAC does have RNA-seq on
  the same 106 patients, but only for the tumours — there are no adjacent
  normal RNA samples, so there's no way to get a tumour vs normal RNA
  comparison out of it. That's why the RNA layer is TCGA, and why this ended up
  being a genuinely cross-cohort analysis.
* **The Vasaikar 2019 supplementary tables wouldn't download.** They have the
  differential results already worked out, which would have saved me a step.
  PMC puts them behind a download challenge, Europe PMC says the article isn't
  open access, and cell.com returns a 403. I calculated the differentials from
  the LinkedOmics matrices instead and checked them against the study's own
  log2FC file.
* **The `cptac` Python package couldn't give me the proteome.** `cptac.Coad()`
  loads transcriptomics and mutations fine, but both proteomics sources
  (`umich` and `bcm`) throw errors in the current version. LinkedOmics had the
  same data as a direct download.
* **Plain mutation frequency was useless on its own.** It put TTN (56.6%),
  MUC16 (36.8%) and OBSCN (27.4%) at or above KRAS (33.0%), which is just a
  gene length effect. I added a CDS-length correction, after which the top of
  the genomic layer is TP53, APC, KRAS, RPL22, ACVR2A and TCF7L2. It's a length
  correction, not a proper driver-gene method, and the report says that.

### One thing about the result

I picked KRAS at the start as the gene I wanted to look at, not as something I
expected to come out on top, and I didn't tune anything to move it. It came out
at rank 821 of 4,851: one of the strongest genomic hits (mutated in 33% of
patients) but flat at RNA and protein, because a G12D mutation changes what the
protein does rather than how much of it there is. The report argues that this
is a limitation of scoring on abundance, not a finding about KRAS.

Packages needed: numpy, pandas, scipy, jupyter (plus biopython for Module 1,
and openpyxl if you want the independent-cohort check).
