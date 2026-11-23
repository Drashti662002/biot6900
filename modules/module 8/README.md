# BIOT 6900 · Module 8: Single-cell and spatial omics (breast cancer)

**Name:** Drashti Bharat Bhanushali
**Target gene:** ERBB2 (HER2)

## Files

| File | What it is |
|---|---|
| `BIOT6900_Module8_Assignment_Starter.ipynb` | **The graded notebook**, Restart & Run All with outputs |
| `Module8_ERBB2.pdf` | **The report** |

## What did not run as designed

- **scVI ran on CPU, not GPU.** The notebook was run on a laptop, so scVI trained for 60 epochs (the notebook's CPU setting) instead of 200. Its scores may improve with full GPU training on Explorer.
- **Census labels needed one correction.** The provided `broad_type()` checks for "t cell" before tumor words, and "malignant cell" contains "t cell", so the Census's 8,526 malignant cells were labeled T. The first TODO cell (1.3) moves them to Epithelial/Tumor before any plot or score.
- **Environment fixes, needed to run on macOS:**
  - The imports cell sets `NUMBA_THREADING_LAYER=workqueue`, because duplicate OpenMP libraries from torch, scikit-learn and conda crashed the kernel.
  - Step 1.2 uses `span=0.5` in `highly_variable_genes`, because two low-depth donors fail the loess fit at the default 0.3.
  - `pynndescent` was pinned to 0.5.13, because 0.7.0 breaks `sc.pp.neighbors` in scanpy 1.11.5.
