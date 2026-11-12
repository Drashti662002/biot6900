# Module 9 — Clinical NLP: signal vs leakage

## Disease
Community-acquired pneumonia (CAP), adult inpatients.

- **Signals (early findings):** confusion, tachypnea, hypoxia (low SpO2 on room air)
- **Consequences (planted leakage):** intubation, vasopressors, comfort_care

## How to re-run
Open `BIOT6900_Module9_Notebook3_Graded.ipynb` and run Kernel → Restart & Run All.

- Needs pandas, scikit-learn and matplotlib (tested with Python 3 from Anaconda, scikit-learn 1.9.0).
- `SEED = 2088061` (student ID 002088061, leading zeros dropped for Python), so results are reproducible.

## Key results

| Feature set | CV AUROC |
|---|---|
| All features | 0.891 ± 0.044 |
| Declared leakage removed | 0.786 ± 0.032 |
| All leakage removed | 0.720 ± 0.060 |

Final leak set: intubation, vasopressors, comfort_care, renal_replacement. I kept lactate_not_clearing as a signal.

## What did not work / limitations
- Every synthetic note says "admitted to the ICU", so this is really an ICU (severe CAP) cohort, not general CAP inpatients.
- Negated notes always use the first synonym, which gives odd phrases like "No evidence of confused".
- The notes have no timestamps, so timing is assumed by the generator. In real charts I could not rely on note order to separate signal from leakage.
