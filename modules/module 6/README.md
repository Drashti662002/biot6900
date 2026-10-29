# BIOT 6900 · Module 6: Polygenic risk scores

**Name:** Drashti Bharat Bhanushali
**Disease (Part 3):** colorectal cancer, carried forward from my Module 2–4 colon adenocarcinoma work
**Score:** PGS Catalog PGS000148 (CRC63, Jeon et al. 2018): 63 SNPs with European-only weights

## Files

| File | What it is |
|---|---|
| `BIOT6900_Module6_Toy_PRS.ipynb` | Part 1, toy model, executed (AUC 0.659, or 0.536 without the APOE-like variant) |
| `BIOT6900_Module6_Part2_RealPRS.ipynb` | Part 2, worked example on the Canvas Alzheimer's data, executed |
| `BIOT6900_Module6_Part3_Graded_Starter.ipynb` | Part 3 starter, unchanged |
| `data/`, `data.zip` | Canvas Alzheimer's data (used by Part 2 only) |
| `part3_crc/BIOT6900_Module6_Part3_Graded_CRC.ipynb` | **The graded notebook**: colorectal cancer, `USE_DEMO = False`, Restart & Run All |
| `part3_crc/module 6_report.docx` | **The report**, with the proof-of-data block |
| `part3_crc/prepare_crc_data.py` | Builds the four input files from the PGS Catalog and 1000 Genomes |
| `part3_crc/data/` | `genotypes.csv`, `variants.csv`, `score.csv`, `meta.csv`, `simulation_log.csv`, `README.md` |
| `part3_crc/figures/` | PRS-by-ancestry and calibration plots |

## Why I built my own data

The Canvas `data.zip` holds only the 22 Alzheimer's SNPs, so a colorectal score has nothing to match against. `prepare_crc_data.py` rebuilds the same four files from the same sources:

- **Score:** the PGS000148 harmonized GRCh37 scoring file.
- **Genotypes:** 1000 Genomes phase 3 (GRCh37), read region by region with bcftools.
- **People:** the same 2,336 people, ancestry labels and ages as the class cohort.
- **Phenotype:** simulated with the professor's recipe from the Canvas README. Every number is unchanged and I did no tuning.

## What didn't work

- **rs202110856** has no record in 1000 Genomes phase 3, so 62 of 63 SNPs are in the genotype file.
- **rs7229639**'s published row has no effect allele. It is left blank, so the notebook's keep/flip/drop rule drops it (used 61, flipped 28, dropped 1).
- **The EBI FTP connection hung once mid-download**, so the script gives each region a 5-minute timeout and retries up to 4 times.

## Rebuilding the data

```bash
conda create -n prs_tools -c conda-forge -c bioconda bcftools   # one time
cd "modules/module 6/part3_crc"
PATH="$(conda info --base)/envs/prs_tools/bin:$PATH" python prepare_crc_data.py
```

The script needs an internet connection and takes about 15 minutes. Cached regions in `part3_crc/data/raw/` are gitignored.
