# BIOT 6900 · Module 4 Assignment: Structure-Based Druggability

**Name:** Drashti Bharat Bhanushali
**Disease:** colon adenocarcinoma (`coad`)
**Targets:** CA2 (UniProt P00918, Module 2 rank 1, PDB 3HS4 chain A) and CR2 (UniProt P20023, Module 2 rank 11, PDB 1LY2 chain A)

`BIOT6900_Module4_Assignment.ipynb` assesses these two targets.

## What didn't work

- **CR2 Step 5 comparison.** 1LY2 covers only CR2 Sushi domains 1–2 (UniProt 20–149), but the AlphaFold model is 1,033 residues long. The provided `global_align` scatters 48 of the 130 pairs outside 20–149, so the printed RMSD (38.45 Å), the "residues that differ most" for CR2, and the CR2 ligand-site mapping are alignment artefacts. I did not edit the provided code. An outside check aligning only AlphaFold residues 20–149 gave 1.00 Å over 121 pairs.
- **fpocket volumes vary slightly between runs.** Scores, pLDDT, PAE and verdicts were identical across two runs; only volumes changed by a few Å³.

## What the notebook does

- **Step 1 · Settings:** disease, the two targets (UniProt accessions), their Module 2 rank and reason, and the PDB entry and chain for each (CA2: 3HS4 chain A; CR2: 1LY2 chain A).
- **Step 2:** downloads each target's AlphaFold model and PAE file.
- **Step 3:** plots per-residue pLDDT and tabulates mean pLDDT, % > 90, % < 50 and stretches of 10+ residues below 50.
- **Step 4:** plots the PAE heatmaps, finds confident domains and computes the mean PAE between domains 1 and 2.
- **Step 5:** pairs AlphaFold and experimental residues by sequence alignment, superimposes on identical, confident (pLDDT ≥ 70) pairs, reports the RMSD (per domain when inter-domain PAE > 10 Å) and plots per-residue deviation.
- **Step 6:** runs fpocket on each AlphaFold model, adds mean pLDDT and PAE among lining residues to each pocket, and labels it reject / suspect / keep. A ligand-site check then compares bound ligands in the experimental structure with the pockets.
- **Step 7:** records a druggability call per target and exports the summary table.

## Environment

- The `biot6900` conda environment (Jupyter kernel `python3` in that environment), with numpy, pandas and matplotlib.
- fpocket, installed in that environment with `conda install -c conda-forge fpocket`.

## How to run

Open the notebook in Jupyter or VS Code with the `biot6900` kernel, then Kernel → Restart & Run All. It needs an internet connection, because it downloads files from the AlphaFold Database and RCSB.

From the command line:

```
conda activate biot6900
jupyter nbconvert --to notebook --execute --inplace BIOT6900_Module4_Assignment.ipynb
```

## Output files

- `structure_assessment_coad.csv`: one row per target (Step 7 summary).
- `Module4_Assignment_Report.pdf`: the report (Steps 8 and 9).
- `structures/`: downloaded AlphaFold models, PAE files, PDB entries and fpocket output. This is not committed (see `.gitignore`).
