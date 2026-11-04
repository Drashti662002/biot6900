# BIOT 6900 · Module 7: Disease knowledge graph (Neo4j)

**Name:** Drashti Bharat Bhanushali
**Disease:** colon adenocarcinoma (COAD), carried forward from Modules 2–6

## Files

| File | What it is |
|---|---|
| `BIOT6900_Module7_Starter.ipynb` | **The graded notebook**, Restart & Run All with outputs |
| `Module7_KnowledgeGraph_Bhanushali.docx` | **The report** |
| `my_knowledge_graph.png` | Graph figure, coloured by node type, with legend |
| `nodes.csv`, `edges.csv` | Node and edge export (every edge has `source` and `confidence`) |
| `prepare_kg_data.py` | Builds `data/targets_coad.csv` from my Module 2, 3 and 6 outputs |
| `data/targets_coad.csv` | The 15 COAD genes with their evidence, read by the notebook |
| `data/gwas_lookup.csv` | Cached GWAS Catalog lookups for the lead SNPs |

## Running it

1. Copy your Aura credentials file to `.env` in this folder (it is git-ignored).
2. Open `BIOT6900_Module7_Starter.ipynb` and run Kernel → Restart & Run All.

`python prepare_kg_data.py` rebuilds `data/targets_coad.csv`. It needs the shared Module 2 and Module 3 data.

The graph has 54 nodes and 61 relationships. Ulcerative colitis is the related disease that drives the repurposing query.
