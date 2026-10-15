# BIOT 6900 coursework — Drashti Bharat Bhanushali

## Repository layout

Each module is a self-contained folder. **Open notebooks from inside their own
module folder** — they use relative paths, so Jupyter's working directory has to
be the module directory.

```
modules/
├── module1/   environment setup, Biopython, database lookups
├── module2/   multi-omics target identification (colon adenocarcinoma)
│   └── data/        -> shared/module2_omics
├── module3/   protein language models: toy, real data, assignment
│   ├── data/        -> shared/module3_pack
│   └── targets_coad.csv  -> ../module2/targets_coad.csv
├── module4/   predicted structures and drug pockets
│   └── structures/  -> shared/module4_structures
└── module5/   variant interpretation from first principles (no data needed)

shared/        one physical copy of every large dataset
docs/          lab guides (PDF)
```

The `data/` and `structures/` entries inside each module are **symlinks** into
`shared/`, so nothing is duplicated on disk and `Path("data")` still resolves
exactly as the notebooks expect. Module 3 reads Module 2's `targets_coad.csv`
through a symlink, so there is only ever one copy of that file.

Running anything:

```bash
conda activate biot6900
cd modules/module3        # or whichever module
jupyter notebook          # then Kernel → Restart & Run All
```

Large inputs are gitignored, with a comment in `.gitignore` saying how to
restore each one. Nothing gitignored is produced by this repo.

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
cd modules/module2
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

## Module 3 — Protein Language Models for Target Triage & Variant Effects

| File | What it is |
|---|---|
| `BIOT6900_Module3_Part1_ToyModel.ipynb` | Lab Part 1, the toy model. Runs offline. |
| `BIOT6900_Module3_Part2_RealData.ipynb` (and `_2`) | Lab Part 2 on ~1,200 real human proteins. Identical apart from one sentence. |
| `BIOT6900_Module3_Assignment.ipynb` | The graded assignment, run on my Module 2 colon adenocarcinoma targets. |
| `targets_coad_w3.csv` | The re-ranked target list handed to Week 4. |
| `data/` | Symlink to the Canvas data pack in `shared/module3_pack/`. |

The data pack zip and the `.npy` embeddings are gitignored. To restore them,
unzip `BIOT6900_Module3_DataPack.zip` into `shared/module3_pack/`.

## Module 4 — Predicted structures and drug pockets

| File | What it is |
|---|---|
| `BIOT6900_Module4_Part1_ToyModel.ipynb` | Toy model for pLDDT, PAE and pockets. |
| `BIOT6900_Module4_Demo_DrugPocket.ipynb` | Acetazolamide (crystal 3HS4) placed in the AlphaFold model of carbonic anhydrase II. |
| `BIOT6900_Module4_Part2_RealData_1.ipynb` | The same checks on p53, APOE and TREM2. |
| `demo_pymol.pml`, `demo_chimerax.cxc` | Scripts that open the aligned demo scene. |

`structures/` is gitignored — the notebooks re-download the AlphaFold models and
PDB entries on first run. Section 2.6 needs `fpocket`
(`conda install -c conda-forge fpocket`).

## Module 5 — Variant interpretation from first principles

`BIOT6900_Module5_Part1_Toy.ipynb` only. No data files, no downloads, no
network: the toy gene is 90 bases written into the notebook.

## Packages

numpy, pandas, scipy, matplotlib, jupyter, and:

* **Module 1** — biopython
* **Module 2** — openpyxl, for the independent-cohort check only
* **Module 3** — scikit-learn
* **Module 4** — py3Dmol, and fpocket for the pocket section
