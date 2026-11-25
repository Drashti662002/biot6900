# BIOT 6900 coursework — Drashti Bharat Bhanushali

Precision-medicine coursework built end to end on one disease: **colon
adenocarcinoma (COAD)**. Every module uses the same disease, so the modules
compose — Module 2's ranked target list feeds Module 3, and the Part II labs
(5–8) are synthesised in the two portfolio deliverables.

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
├── module5/   variant interpretation and pharmacogenomics
├── module 6/  polygenic risk scores (part3_crc = the graded colorectal run)
├── module 7/  Neo4j knowledge graph of COAD targets, drugs and variants
├── module 8/  single-cell and spatial transcriptomics (ERBB2)
├── module 9/  clinical NLP and label leakage
└── module 10/ medical imaging and confounder detection

shared/          one physical copy of every large dataset
docs/            lab guides (PDF)
deliverable_1/   Portfolio Deliverable 1 — Part I target triage
deliverable_2/   Portfolio Deliverable 2 — Part II COAD case study
deliverables/    the individual Lab 5-8 reports
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

## Module 5 — Variant interpretation and pharmacogenomics

`BIOT6900_Module5_Part1_Toy.ipynb` needs nothing: the toy gene is 90 bases
written into the notebook. Part 2 and the graded Part 3 pull real variants from
MyVariant.info and score them with CADD, REVEL and AlphaMissense. The graded run
calibrates CADD across five pharmacogenes (DPYD, CYP2C19, CYP2C9, TPMT, VKORC1)
and reports a positive likelihood ratio per score band.

`modules/module5/data/*.tsv` are the cached MyVariant.info fetches, one per
gene. They are committed because Deliverable 2 quotes per-variant scores out of
them; delete one and re-run the notebook to refetch it.

## Module 6 — Polygenic risk scores

The graded colorectal run is `part3_crc/BIOT6900_Module6_Part3_Graded_CRC.ipynb`:
PGS000148 applied to 1000 Genomes genotypes across five ancestries, with
discrimination, calibration and tag-causal r² per ancestry.

**The phenotype in this lab is simulated.** It is generated from the score, so
the AUCs measure the simulation, not clinical performance. Every figure and
claim derived from it is flagged as simulated.

```bash
cd "modules/module 6/part3_crc" && python prepare_crc_data.py   # needs bcftools
```

## Module 7 — Knowledge graph

`BIOT6900_Module7_Starter.ipynb` builds a 54-node, 61-edge COAD graph in Neo4j
AuraDB, then runs PageRank, greedy-modularity communities and multi-hop Cypher.
Every edge carries a source and a confidence. `nodes.csv`, `edges.csv`,
`data/targets_coad.csv` and `data/gwas_lookup.csv` are committed, so the
analysis is reproducible without a database.

### Neo4j credentials

Credentials are **never committed**. Put them in `modules/module 7/.env`, which
is gitignored along with `Neo4j-*.txt` and `*.credentials`:

```bash
NEO4J_URI=neo4j+s://YOUR-INSTANCE.databases.neo4j.io
NEO4J_USERNAME=YOUR-INSTANCE-ID     # the instance ID, not "neo4j"
NEO4J_PASSWORD=YOUR-PASSWORD
```

The notebooks read these with `os.environ.get`, and the values checked into the
notebooks are placeholders.

## Module 8 — Single-cell and spatial transcriptomics

`BIOT6900_Module8_Assignment_Starter.ipynb`: ERBB2 across cell types in a CZ
CELLxGENE Census breast-cancer subset, integrated with Harmony and scVI and
compared by silhouette score, then a 10x Visium section analysed with Squidpy
(Moran's I, neighbourhood enrichment).

## Module 9 — Clinical NLP

Label leakage in community-acquired pneumonia notes: the phrase that predicts
the label is the phrase a clinician writes *because* of the label.

## Module 10 — Imaging and confounder detection

`BIOT6900_Module10_Confounder_Experiment.ipynb` trains a MONAI DenseNet121 on
PathMNIST tumour vs normal, then injects a corner artifact into one class only
and tests whether Grad-CAM, occlusion and a clean control model can detect the
shortcut. CPU only. The raw `.npz` is never modified — arrays are copied in
memory.

## Portfolio deliverables

| Folder | What it is |
|---|---|
| `deliverable_1/` | Part I target triage: the 15-gene evidence table and its write-up. |
| `deliverable_2/` | Part II COAD case study synthesising Labs 5–8 and Deliverable 1. `Portfolio_deliverable_Drashti.pdf` is the submitted document; `inventory.md` and `provenance.md` trace every number to the file or notebook cell it came from. |

Deliverable 2 performs no new analysis: `build_evidence_table.py` and
`make_figures.py` only read saved lab outputs and re-plot them.

## Data sources

All open access, no data use agreement.

| Source | Used in |
|---|---|
| TCGA-COAD (UCSC Xena GDC hub) | Module 2 RNA layer |
| CPTAC-2 Prospective Colon Cancer (LinkedOmics) | Module 2 protein and mutation layers |
| ClinVar and gnomAD, via MyVariant.info | Module 5 |
| PGS Catalog (PGS000148) | Module 6 |
| 1000 Genomes Phase 3 | Module 6 genotypes |
| GWAS Catalog | Modules 1, 6, 7 |
| ChEMBL, Reactome, HPO | Module 7 edges |
| CZ CELLxGENE Census | Module 8 single-cell |
| 10x Genomics Visium (human breast cancer, block A section 1) | Module 8 spatial |
| MedMNIST PathMNIST | Module 10 |

## Packages

numpy, pandas, scipy, matplotlib, jupyter, and:

* **Module 1** — biopython
* **Module 2** — openpyxl, for the independent-cohort check only
* **Module 3** — scikit-learn
* **Module 4** — py3Dmol, and fpocket for the pocket section
