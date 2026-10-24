# BIOT 6900 · Module 5 Assignment: Variant Interpretation

**Name:** Drashti Bharat Bhanushali
**Gene:** CYP2C19
**AlphaGenome Atlas:** not run — no API key, so `SKIP_ATLAS = True` and the Atlas part of the report is answered from reasoning

`BIOT6900_Module5_Part3_Assignment.ipynb` runs the full workflow on CYP2C19: fetch, coverage audit, frequency filter, calibration, VUS classification and discordance.

## What didn't work

- **The assignment notebook would not execute as distributed.** Every markdown cell had `cell_type: "md"`, which is not valid under the nbformat schema, and no cell carried an `id`. `jupyter nbconvert --execute` refused the file. I corrected both (`md` → `markdown`, added ids) without changing any cell content, and the notebook then ran top to bottom.
- **CYP2C19 cannot be calibrated on its own.** After filtering it has **0 P/LP** and 5 B/LB variants with a CADD score, because ClinVar files this gene's functional alleles as *drug response* rather than pathogenic. There are no positives to count, so every likelihood ratio would be 0 or undefined. I took the documented route and pooled with the other four lab genes (CYP2C9, VKORC1, DPYD, TPMT) using the same fetch and the same 0.5 filter, giving 167 labelled variants. 39 of the 49 positives come from DPYD, so the result is largely a DPYD calibration transferred to CYP2C19. This is disclosed in the notebook output and discussed in section 1 of the report.
- **No AlphaGenome key**, so Step 7 was skipped. The cost is specific: the 25 non-missense variants in this gene are exactly what the Atlas is for.
- **`rs4244285` (\*2) comes back as two rows**, `chr10:g.94781859G>A` and `G>T`, two different alternate bases at the same position. Not an error, but the anchor check reports both.

## What the notebook does

- **Step 0 · Setup:** API fields, parsing helpers, `fetch_gene()`, `numeric()`, `clinvar_label()` and `calibrate()`, as provided. `MY_GENE = "CYP2C19"`.
- **Q1 · Fetch and quality-check:** pulls every ClinVar variant for the gene from MyVariant.info, caches to `data/CYP2C19.tsv`, then asserts on three checks — variant count (80), both anchor alleles present (\*2, \*17), and the score columns non-empty (CADD 85%, REVEL and AlphaMissense 67.5%, gnomAD 68.8%).
- **Q2 · Coverage audit:** scores per tool broken down by consequence class, the overall fraction each tool could score, and which classes have zero REVEL scores (8 classes, 25 variants).
- **Q3 · Frequency filter:** `MAX_AF = 0.5`, chosen for a pharmacogene rather than a rare-disease gene. Reports what it removes (2 variants, both benign), confirms no functional allele is lost, and contrasts against a 0.001 filter (19 removed, 14 of them drug-response or anchor alleles). Variants with no gnomAD frequency are kept.
- **Q4 · Calibration:** calibrates CADD over bands `[0, 10, 20, 25, 30, 50]`, pooling as described above. Reports the likelihood ratio and ACMG evidence strength per band, the number of labelled variants behind the whole thing, which bands rest on fewer than 10, and a Haldane-corrected LR for the band whose LR is infinite because a cell is empty.
- **Q5 · VUS classification:** maps each VUS to its band's evidence, but only where that band rests on at least 10 labelled variants. Exports `pgx_classified_CYP2C19.csv`.
- **Q6 · Atlas and discordance:** the Atlas cell is skipped; the discordance cell applies the Part 2 rules and flags 21 variants.

## Environment

- The `biot6900` conda environment (Jupyter kernel `python3` in that environment), with numpy, pandas, matplotlib and requests.
- An internet connection is needed on the first run. Everything then caches to `data/`.

## How to run

Open the notebook in Jupyter or VS Code with the `biot6900` kernel, then Kernel → Restart & Run All.

From the command line:

```
conda activate biot6900
jupyter nbconvert --to notebook --execute --inplace BIOT6900_Module5_Part3_Assignment.ipynb
```

If a download stalls partway, delete `data/CYP2C19.tsv` before retrying. A partial cache file is loaded without complaint and is worse than none.

## Output files

- `pgx_classified_CYP2C19.csv`: the 33 VUS with their calibrated evidence and provisional call. Of these, 8 lean pathogenic, 16 lean benign and 9 do not move (3 in a band too small to apply, 6 large deletions with no score).
- `Module5_Assignment_Report_final.docx`: the report.
- `data/`: cached MyVariant.info fetches, one TSV per gene. Not committed — the notebook re-downloads them on first run.
