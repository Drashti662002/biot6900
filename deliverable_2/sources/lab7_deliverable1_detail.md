# Evidence inventory: Lab 7 (COAD knowledge graph) + Portfolio Deliverable 1

Read-only extraction, 24 Nov 2026. No Neo4j connection, no re-execution, no `.env` contents read.

Key source paths (all absolute):
- `/Users/drashtibhanushali/biot6900/modules/module 7/nodes.csv`, `edges.csv`, `prepare_kg_data.py`, `README.md`, `my_knowledge_graph.png`, `Module7_report.pdf`, `Module7_report.docx`, `BIOT6900_Module7_Starter.ipynb`
- `/Users/drashtibhanushali/biot6900/modules/module 7/data/targets_coad.csv`, `data/gwas_lookup.csv`
- `/Users/drashtibhanushali/biot6900/deliverable_1/deliverable1.md`, `Deliverable1_Bhanushali.docx/.pdf`, `evidence_table.csv`, `inventory.md`, `provenance.md`
- `/Users/drashtibhanushali/biot6900/modules/module3/targets_coad_w3.csv`, `/Users/drashtibhanushali/biot6900/modules/module2/targets_coad.csv`

**Important notebook caveat.** The working-tree copy of `BIOT6900_Module7_Starter.ipynb` (modified, uncommitted) has been reverted to a *partly un-run* state: cell 3 fails ("Not connected yet — Failed to DNS resolve YOUR-INSTANCE.databases.neo4j.io"), cell 7 errors (`ServiceUnavailable`), and cells 10–18 have **no outputs** and no execution counts. It also lacks the `CONFIDENCE` dict, so `edges.csv` could not have been produced by it. All executed query / PageRank / community / figure outputs quoted below come from the **committed** version at `git show HEAD:"modules/module 7/BIOT6900_Module7_Starter.ipynb"` (commit `bf72daa`), which is fully run (exec 1–10, "Connected to Neo4j", node counts 54 / 61 returned by the live database).

---

## 1. Node and edge counts by type, with source and confidence fields

`nodes.csv` columns: `name, kind, ensembl, gwas_p, rna_lfc, prot_lfc, druggability, status` (54 data rows).
`edges.csv` columns: `src, rel, dst, confidence, pvalue, source` (61 data rows). Every edge has both `source` and `confidence`; `pvalue` only on GWAS gene–disease edges.

### Nodes by kind (counted from nodes.csv)

| kind | n | members |
|---|---|---|
| gene | 18 | CA2, DES, CR2, MYH11, KRAS, BRAF, EGFR, VEGFA, TYMS, TOP1, PTGS2, MYC, SMAD7, BMP4, TCF7L2 (15 COAD) + IL23R, JAK2, IL12B (UC) |
| disease | 2 | Colon adenocarcinoma, Ulcerative colitis |
| phenotype | 5 | Hematochezia, Abdominal pain, Iron deficiency anemia, Weight loss, Constipation |
| pathway | 9 | WNT/beta-catenin; TGF-beta/BMP; EGFR-RAS-MAPK; VEGF signaling / angiogenesis; Nucleotide synthesis & DNA replication; Prostaglandin synthesis; Carbon dioxide hydration; Complement cascade; Smooth muscle contraction |
| drug | 13 | Cetuximab, Panitumumab, Bevacizumab, Encorafenib, Sotorasib, Adagrasib, Fluorouracil, Capecitabine, Irinotecan, Aspirin, Acetazolamide, Ustekinumab, Tofacitinib |
| variant | 7 | rs6983267, rs4939827, rs4444235, rs11196172, KRAS G12D, KRAS G12C, BRAF V600E |
| **total** | **54** | matches report ("54 nodes and 61 relationships") and the Neo4j count in cell 7 |

Node `status` is populated only on drug nodes: 9 × "approved (CRC)", 2 × "approved (other indication)" (Aspirin, Acetazolamide), 2 × "approved (UC)" (Ustekinumab, Tofacitinib).

### Edges by relationship × source × confidence (counted from edges.csv)

| rel | source | confidence | n | pvalue present |
|---|---|---|---|---|
| ASSOCIATED_WITH | multi-omics (TCGA/CPTAC) | medium | 4 (CA2, DES, CR2, MYH11 → COAD) | no |
| ASSOCIATED_WITH | somatic mutation (CPTAC) | high | 2 (KRAS, BRAF → COAD) | no |
| ASSOCIATED_WITH | pharmacology | low | 5 (EGFR, VEGFA, TYMS, TOP1, PTGS2 → COAD) | no |
| ASSOCIATED_WITH | GWAS | high | 4 (MYC 2e-122, SMAD7 3e-30, BMP4 8e-10, TCF7L2 7e-20 → COAD) | yes |
| ASSOCIATED_WITH | GWAS Catalog | high | 3 (IL23R 5e-28, JAK2 2e-25, IL12B 1e-21 → Ulcerative colitis) | yes |
| HAS_PHENOTYPE | HPO | low | 8 (5 from COAD, 3 from UC) | no |
| PARTICIPATES_IN | Reactome | low | 15 | no |
| TARGETS | ChEMBL | high | 13 | no |
| OCCURS_IN | GWAS Catalog (Module 6 PRS) | low | 4 | no |
| OCCURS_IN | CPTAC somatic (Module 2) | high | 3 | no |
| **total** | | | **61** | |

Confidence totals: high 25, medium 4, low 32.

Confidence is assigned by the `CONFIDENCE` dict in cell 5 of the **committed** notebook (not present in the working copy), keyed on the edge's `source`:
```
GWAS: high                       # genome-wide significant (p < 5e-8) CRC locus
GWAS Catalog: high               # genome-wide significant UC locus
somatic mutation (CPTAC): high   # recurrent driver mutation in CPTAC tumours
multi-omics (TCGA/CPTAC): medium # concordant RNA + protein change; tissue composition confounds it
pharmacology: low                # in the graph because a drug exists (circular)
HPO: low                         # hand-entered, nonspecific presenting symptoms
Reactome: low                    # my own simplified pathway groupings
ChEMBL: high                     # established drug mechanism of action
CPTAC somatic (Module 2): high   # coding hotspot inside the gene
GWAS Catalog (Module 6 PRS): low # non-coding lead SNP assigned to a candidate gene
```
Report Table 1 (Module7_report, §1) gives the same per-relationship counts: ASSOCIATED_WITH 18, HAS_PHENOTYPE 8, PARTICIPATES_IN 15, TARGETS 13, OCCURS_IN 7.

---

## 2. PageRank table (full, with values)

Source: committed notebook cell 14 output (`nx.pagerank` on the directed graph; `rounded to 4 dp`). Columns as printed.

| gene | pagerank | gwas_p | druggability |
|---|---|---|---|
| KRAS | 0.0377 | NaN | 0.498057 |
| BRAF | 0.0231 | NaN | 0.690298 |
| EGFR | 0.0231 | NaN | 0.708745 |
| TYMS | 0.0231 | NaN | 0.093201 |
| CA2 | 0.0158 | NaN | 0.769118 |
| MYC | 0.0158 | 2.000000e-122 | 0.028789 |
| JAK2 | 0.0158 | 2.000000e-25 | 0.673000 |
| TCF7L2 | 0.0158 | 7.000000e-20 | 0.038543 |
| BMP4 | 0.0158 | 8.000000e-10 | 0.079581 |
| SMAD7 | 0.0158 | 3.000000e-30 | 0.022749 |
| TOP1 | 0.0158 | NaN | 0.091370 |
| PTGS2 | 0.0158 | NaN | 0.463987 |
| VEGFA | 0.0158 | NaN | 0.015645 |
| IL12B | 0.0158 | 1.000000e-21 | 0.716000 |
| DES | 0.0086 | NaN | 0.001530 |
| MYH11 | 0.0086 | NaN | 0.088702 |
| IL23R | 0.0086 | 5.000000e-28 | 0.463000 |
| CR2 | 0.0086 | NaN | 0.640310 |

Report §2 Table 2 compresses this into ranks: 1 KRAS 0.0377; 2= BRAF / EGFR / TYMS 0.0231; 5= CA2, MYC, SMAD7, TCF7L2, BMP4, TOP1, PTGS2, VEGFA 0.0158 (druggability range given as 0.02–0.77); last DES, MYH11, CR2 0.0086 (0.00–0.64). The report's Table 2 covers only the 15 COAD genes — JAK2 (0.0158) and IL23R (0.0086) from the output are omitted from the report table.

Report's own reading: only drug and variant edges point *into* genes, so PageRank = "how many drugs and variants I attached"; floor value 0.0086.

---

## 3. Community detection (membership per community)

Source: committed notebook cell 16, `greedy_modularity_communities(G.to_undirected())` → **10 communities**:

| # | members |
|---|---|
| 0 | Abdominal pain, Hematochezia, IL12B, IL23R, Iron deficiency anemia, JAK2, Tofacitinib, Ulcerative colitis, Ustekinumab |
| **1** | **CR2, Colon adenocarcinoma, Complement cascade, Constipation, DES, MYH11, Smooth muscle contraction, Weight loss** |
| 2 | BRAF, BRAF V600E, Cetuximab, EGFR, EGFR-RAS-MAPK signaling, Encorafenib, Panitumumab |
| 3 | Capecitabine, Fluorouracil, Irinotecan, Nucleotide synthesis & DNA replication, TOP1, TYMS |
| 4 | Adagrasib, KRAS, KRAS G12C, KRAS G12D, Sotorasib |
| 5 | MYC, TCF7L2, WNT/beta-catenin signaling, rs11196172, rs6983267 |
| 6 | BMP4, SMAD7, TGF-beta/BMP signaling, rs4444235, rs4939827 |
| 7 | Bevacizumab, VEGF signaling / angiogenesis, VEGFA |
| 8 | Aspirin, PTGS2, Prostaglandin synthesis |
| 9 | Acetazolamide, CA2, Carbon dioxide hydration |

**Community containing DES, MYH11, CR2 = community 1.** It also contains the COAD disease node, the "Smooth muscle contraction" and "Complement cascade" pathways, and two phenotypes (Constipation, Weight loss). **CA2 is NOT in it** — CA2 sits alone in community 9 with Acetazolamide and Carbon dioxide hydration.

Report interpretation (§2, final paragraph, verbatim sense): "Community detection put DES, MYH11 and CR2 in one cluster with the COAD node, the 'Smooth muscle contraction' and 'Complement cascade' pathways and two symptoms. I read this cluster as **tissue composition, not tumour mechanism**. DES and MYH11 are smooth-muscle proteins and CR2 is a B-cell receptor. Normal colon samples contain muscle and lymphoid tissue, while tumour samples are mostly epithelium, so these genes look 'lost in tumour' partly because the normal samples contain more muscle. Read this way, the cluster is a warning about my own Module 2 ranking." §4 adds: multi-omics edges "may be measuring tissue composition, and the RNA and protein layers come from different patients… I would expect DES and MYH11 to lose their edges."

---

## 4. The 2-hop drug result (query + results)

Query (committed notebook cell 10):
```cypher
MATCH (d:Drug)-[:TARGETS]->(g:Gene)-[:ASSOCIATED_WITH]->(dis:Disease)
WHERE dis.name = $disease
RETURN d.name AS drug, g.name AS gene, d.status AS status
ORDER BY drug
```
Result (11 rows, exec 5):

| drug | gene | status |
|---|---|---|
| Acetazolamide | CA2 | approved (other indication) |
| Adagrasib | KRAS | approved (CRC) |
| Aspirin | PTGS2 | approved (other indication) |
| Bevacizumab | VEGFA | approved (CRC) |
| Capecitabine | TYMS | approved (CRC) |
| Cetuximab | EGFR | approved (CRC) |
| Encorafenib | BRAF | approved (CRC) |
| Fluorouracil | TYMS | approved (CRC) |
| Irinotecan | TOP1 | approved (CRC) |
| Panitumumab | EGFR | approved (CRC) |
| Sotorasib | KRAS | approved (CRC) |

5-hop repurposing query (cell 11) returned 2 rows: **Tofacitinib** (JAK2, via Ulcerative colitis, shared symptoms Hematochezia / Abdominal pain / Iron deficiency anemia) and **Ustekinumab** (IL12B, same route).

---

## 5. The IL12B path (nodes, edges, source + confidence) and the IL-23 hypothesis

Path query (cell 12, committed version returns `type(r) + ' [' + r.source + ', ' + r.confidence + ']'`), output:

```
['Ustekinumab', 'IL12B', 'Ulcerative colitis', 'Iron deficiency anemia', 'Colon adenocarcinoma']
['TARGETS [ChEMBL, high]',
 'ASSOCIATED_WITH [GWAS Catalog, high]',
 'HAS_PHENOTYPE [HPO, low]',
 'HAS_PHENOTYPE [HPO, low]']
```

Edge-by-edge, from `edges.csv`:

| # | edge | source | confidence | pvalue |
|---|---|---|---|---|
| 1 | Ustekinumab –TARGETS→ IL12B | ChEMBL | high | — |
| 2 | IL12B –ASSOCIATED_WITH→ Ulcerative colitis | GWAS Catalog | high | 1e-21 |
| 3 | Ulcerative colitis –HAS_PHENOTYPE→ Iron deficiency anemia | HPO | low | — |
| 4 | Colon adenocarcinoma –HAS_PHENOTYPE→ Iron deficiency anemia | HPO | low | — |

Node properties: IL12B gene, ENSG00000113302, gwas_p 1e-21, druggability 0.716, no rna_lfc/prot_lfc (COAD omics not applicable to the UC genes — `prepare_kg_data.py` sets them NaN; the 3 UC genes are hand-entered in the notebook's `RELATED` dict, not in `targets_coad.csv`). Ustekinumab status "approved (UC)". Lead SNP in `gwas_lookup.csv`: rs6871626, UC, p 1e-21, PMID 21297633, GCST000964, Anderson CA, mapped gene IL12B-AS1.

**Report §3 hypothesis (p40 vs p19 discussion present):** "IL12B codes for p40, the subunit shared by IL-12 and IL-23. IL-23 drives gut inflammation, and UC itself raises colorectal cancer risk. My hypothesis is that blocking IL-23 reduces the growth of colorectal tumours that start on an inflammatory background. One complication the graph cannot see: ustekinumab blocks p40, so it also blocks IL-12, which helps the immune system attack tumours. A drug blocking only the IL-23-specific p19 subunit, such as **mirikizumab**, would be a cleaner test."
Test plan: (1) compare an IL-23/IL-17 signature (IL23A, IL12B, IL23R, IL17A) tumour vs normal in TCGA-COAD ("data I already have"); (2) AOM/DSS colitis-associated cancer mouse model, three arms (anti-p40, anti-p19, control antibody), count tumours.
Report also calls the path weak: the two diseases are linked only by three general symptoms, "so taken literally it would suggest a UC drug for any disease that makes people bleed. I treat the path as a pointer and take the hypothesis from the gene at the end of it." Second candidate tofacitinib via JAK2 by the same route.

Note: the report's inline path string renders the edges as `[ChEMBL]`, `[GWAS Catalog, p = 1 × 10⁻²¹]`, `[HPO]`, `[HPO]` — i.e. it drops the confidence tokens the notebook output actually printed and substitutes the p-value.

---

## 6. Pharmacology edges (circularity) and hand-entered / low-confidence phenotype edges

**Pharmacology edges (5, all ASSOCIATED_WITH → Colon adenocarcinoma, confidence `low`):** EGFR, VEGFA, TYMS, TOP1, PTGS2. The notebook's own comment on `"pharmacology": "low"` reads "in the graph because a drug exists (circular)".
Report §4, first limitation block, verbatim sense: "**Pharmacology gene-to-disease edges (5). The weakest, and circular.** I linked EGFR, VEGFA, TYMS, TOP1 and PTGS2 to COAD because a COAD drug targets them, so when the 2-hop query reports that COAD drugs act on COAD genes, part of that answer is something I put in. These edges also drive the ranking in Section 2. With Open Targets association scores instead, I would expect TYMS and TOP1 to drop out." So yes — explicitly called circular.

**Phenotype edges (8, all HPO, all `low`), explicitly hand-entered:** §4 "Phenotype edges (8). **I typed these in by hand from the HPO term list**, and all five COAD phenotypes are very general. My whole repurposing result passes through them. With proper HPO annotations, or more specific phenotypes such as chronic mucosal inflammation, UC would probably still connect but the symptom-only paths would disappear." Rows: COAD → Hematochezia / Abdominal pain / Iron deficiency anemia / Weight loss / Constipation; UC → Hematochezia / Abdominal pain / Iron deficiency anemia.

**Other weak classes flagged in §4:** GWAS variant→gene edges (4, `low`) — "OCCURS_IN is simply wrong for at least two. rs6983267 sits in an enhancer the Catalog maps to CASC8/CCAT2, and rs4444235 maps near BMP4, not in it… the edge should be called REGULATES." Pathway edges (15, labelled Reactome, `low`) — "they are my own simplified groupings, and four of my nine pathways hold a single gene." Multi-omics edges (4, `medium`) — real but possibly tissue composition, cross-cohort. Drug-target + somatic-variant edges (20, `high`) — strongest, but "the catch is direction: acetazolamide inhibits CA2, but CA2 is lost in tumour, so the edge is true and still not useful for therapy." Global: "every edge counts the same. PageRank treats a pathway edge I typed in by hand exactly like a GWAS hit at 10⁻¹²². A weighted PageRank using −log₁₀ p for genetic edges would fix this, and I would expect MYC and SMAD7 to rise and TYMS to fall."

---

## 7. Graph visualisation source

`my_knowledge_graph.png` (409,609 bytes, 3 Nov 20:30, untracked) was produced by **cell 18 of the committed notebook**, which saved it (`plt.savefig("my_knowledge_graph.png", dpi=150)`) in the same cell that wrote `nodes.csv` and `edges.csv` ("saved my_knowledge_graph.png, nodes.csv, edges.csv"; the stored output shows `<Figure size 1300x1100 with 1 Axes>` + a PNG).

Committed cell-18 layout (not in the working copy, which has the older `spring_layout(seed=7, k=0.7)`, figsize 13×8, no legend):
- `KIND_COLOR = {gene:#D41B2C, disease:#3E6F99, phenotype:#2FA3D6, pathway:#00A58C, drug:#9A5FC0, variant:#FF854F, patient:#8A94A3}`
- Layered layout: `LAYER = {drug:0, variant:0, gene:1, disease:2, pathway:2, phenotype:3}`, 10 barycentre sweeps to reduce crossings; `pos = (L, 1-(i+0.5)/len(col))`; figsize 13×11; node_size 650; edges `#C2CCD6`; labels font 7; legend of node types below the axes; patients excluded.

Rendered figure (confirmed by viewing the PNG): title "Colon adenocarcinoma knowledge graph", four columns left→right = drugs+variants, genes, diseases+pathways, phenotypes, with a 6-entry legend (gene, disease, phenotype, pathway, drug, variant).

**Is node/edge data enough to redraw it?** Yes for topology and colouring — `nodes.csv` has `name` + `kind` for all 54 nodes and `edges.csv` all 61 edges, which is everything the layout code consumes (positions are deterministic given the code: fixed LAYER + barycentre sweeps, no RNG; the older spring layout used `seed=7`). The only thing not recoverable from the CSVs alone is the committed layout code itself, which exists in git (`bf72daa`) but **not** in the working-tree notebook.

---

## 8. Drugs attached to each gene (TARGETS edges, all ChEMBL / high)

| gene | drugs |
|---|---|
| TYMS | Fluorouracil, Capecitabine (notebook comment: "oral 5-FU prodrug") |
| TOP1 | Irinotecan |
| EGFR | Cetuximab, Panitumumab |
| KRAS | Sotorasib ("KRAS G12C, with panitumumab"), Adagrasib ("KRAS G12C, with cetuximab") |
| BRAF | Encorafenib ("BRAF V600E, with cetuximab") |
| VEGFA | Bevacizumab |
| PTGS2 | Aspirin ("CRC chemoprevention evidence", status approved (other indication)) |
| CA2 | Acetazolamide ("Module 4 drug-site ligand", approved (other indication)) |
| IL12B | Ustekinumab ("anti-IL-12/23 p40", approved (UC)) |
| JAK2 | Tofacitinib ("pan-JAK inhibitor", approved (UC)) |
| all others (CR2, DES, MYH11, MYC, SMAD7, BMP4, TCF7L2, IL23R) | none |

---

## 9. rs6983267, SMAD7 / rs7229639 and the GWAS lookup rows

**rs6983267: yes, present as a variant node.** Edge `rs6983267 –OCCURS_IN→ MYC`, source "GWAS Catalog (Module 6 PRS)", confidence **low**, no pvalue on the edge. The p-value lives on the *gene* node and the gene–disease edge: `MYC` node `gwas_p = 1.9999999999999998e-122`; edge `MYC –ASSOCIATED_WITH→ Colon adenocarcinoma`, source GWAS, confidence high, `pvalue = 1.9999999999999998e-122` (i.e. the expected **p = 2e-122**; the float noise is a CSV round-trip of 2e-122 — `data/targets_coad.csv` and `data/gwas_lookup.csv` both store the clean `2e-122`). It is in community 5 with MYC, TCF7L2, rs11196172 and WNT/beta-catenin signaling.

**SMAD7: yes**, as a gene node (ENSG00000101665, gwas_p 3e-30, rna_lfc −0.774, druggability 0.022749) with `SMAD7 –ASSOCIATED_WITH→ COAD` (GWAS, high, pvalue 3e-30) and `rs4939827 –OCCURS_IN→ SMAD7` (GWAS Catalog (Module 6 PRS), low).

**rs7229639: `[MISSING: rs7229639 does not appear anywhere in Module 7 — not in nodes.csv, edges.csv, gwas_lookup.csv, targets_coad.csv, prepare_kg_data.py or the report]`.** The only SMAD7 variant in the graph is rs4939827. rs7229639 exists only in Module 6: `/Users/drashtibhanushali/biot6900/modules/module 6/part3_crc/data/variants.csv` line 54 (`rs7229639,A,G`), `score.csv` line 54 (`rs7229639,,G,0.023` — empty effect allele), `simulation_log.csv` line 54 (`rs7229639,0.023,18:46452637,…`, chr18:46,452,637 = the SMAD7 locus), and the Module 6 README: "**rs7229639**'s published row has no effect allele. It is left blank, so the notebook's keep/flip/drop rule drops it (used 61, flipped 28, dropped 1)." So it was dropped from the Module 6 PRS and never carried into the KG.

### `data/gwas_lookup.csv` — all 7 rows verbatim

| gene | gwas_snp | gwas_trait | gwas_p | gwas_pmid | gwas_study | gwas_author | gwas_mapped_genes |
|---|---|---|---|---|---|---|---|
| MYC | rs6983267 | colorectal cancer | 2e-122 | 36539618 | GCST90129505 | Fernandez-Rozadilla C | CASC8;CCAT2;PCAT1;POU5F1B |
| SMAD7 | rs4939827 | colorectal cancer | 3.0000000000000003e-30 | 29917119 | GCST006131 | Schmit SL | SMAD7 |
| BMP4 | rs4444235 | colorectal cancer | 8e-10 | 19011631 | GCST000270 | Houlston RS | DDHD1-DT;MIR5580 |
| TCF7L2 | rs11196172 | colorectal cancer | 7e-20 | 36539618 | GCST90129505 | Fernandez-Rozadilla C | TCF7L2 |
| IL23R | rs11209026 | ulcerative colitis | 5e-28 | 21297633 | GCST000964 | Anderson CA | IL23R |
| JAK2 | rs10758669 | ulcerative colitis | 2e-25 | 21297633 | GCST000964 | Anderson CA | HNRNPA1P41;JAK2 |
| IL12B | rs6871626 | ulcerative colitis | 1e-21 | 21297633 | GCST000964 | Anderson CA | IL12B-AS1 |

Lookup method (`prepare_kg_data.py::gwas_best`): GWAS Catalog REST v2 `associations?rs_id=…`, keeping only associations where **every** `efo_trait` matches the trait regex (so pooled cross-disease studies are excluded), then taking the minimum p-value. Cached in `gwas_lookup.csv`; the cache is reused if present. The report notes the mapped-gene mismatch explicitly (rs6983267 → CASC8/CCAT2 vs MYC attached; rs4444235 → DDHD1-DT/MIR5580, "maps near BMP4, not in it").

---

## 10. Deliverable 1: CA2/CR2 shortlist, the 15-gene triage list, ranking logic, direction problem

### The shortlist (2 targets, not 3–5)
`deliverable1.md` §5 Table 1 and the final `Deliverable1_Bhanushali.docx` (which expands the table with the 13 non-shortlisted genes):

| Rank | Protein | Layers | Concordance (RNA / protein log2FC) | Prior | Mean pLDDT | Pocket pLDDT | Pocket score | Rationale |
|---|---|---|---|---|---|---|---|---|
| 1 | CA2 | 2 of 3 | concordant, down (−5.20 / −2.16) | p_sm 0.770 | 97.4 | 97.9 | 0.298 (pocket 2, acetazolamide site) | confident structure + known-ligand pocket; small-molecule route |
| 2 | CR2 | 2 of 3 | concordant, down (−1.97 / −1.41) | p_ab 0.663 | 74.9 | 65.2 | 0.898 (**rejected**, pocket 50) | surface receptor; antibody route; no confident small-molecule pocket |
| – | TSPAN7 | 2 of 3 | down / down | p_sm 0.754 | n/a | n/a | n/a | not shortlisted: good prior, not structurally assessed |
| – | BCHE | 2 of 3 | down / down | p_sm 0.736 | n/a | n/a | n/a | same |
| – | CFL2 | 2 of 3 | down / down | p_sm 0.525 | n/a | n/a | n/a | moderate prior, not assessed |
| – | DES, PRPH | 2 of 3 | down / down | 0.001 (both) | n/a | n/a | n/a | lowest priors, intermediate filaments, no catalytic site |
| – | other 8 | 2 of 3 | down / down | see Lab 3 | n/a | n/a | n/a | not structurally assessed |

Supporting numbers: CA2 RNA q 2.7e-22, protein q 9.5e-16, mutated 5/106 q 0.50, AF-vs-3HS4 RMSD 0.37 Å over 257 pairs, 98.5% residues > 90. CR2 RNA q 1.7e-13, protein q 2.4e-6, mutated 10/106 q 0.86, 1.3% > 90, 6.3% < 50, 8 confident domains, inter-domain PAE 26.3 Å, best *confident* pocket = #27 score 0.005 pLDDT 75.5, Sushi domain 1 vs 1LY2 2.20 Å over 80 residues, whole-chain 38.45 Å called an alignment artefact.

### The 15-gene triage list (ranks and scores)
Lab 2 rank (`modules/module2/targets_coad.csv`) + Lab 3 re-rank and priors (`modules/module3/targets_coad_w3.csv`), as carried into `deliverable_1/evidence_table.csv`:

| Lab2 rank | gene | UniProt | RNA log2FC (q) | protein log2FC (q) | mutation | Lab3 rank | week2_score | p_sm | p_ab | best prior | final |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | CA2 | P00918 | −5.20 (2.7e-22) | −2.16 (9.5e-16) | 4.7% (5/106) q 0.50 | 1 | 0.9917 | 0.7698 | 0.0085 | 0.770 sm | 1.0000 |
| 2 | DES | P17661 | −4.21 (3.2e-15) | −2.79 (5.2e-15) | 5.7% (6/106) q 0.77 | 6 | 0.9887 | 0.0010 | 0.0014 | 0.001 | 0.6133 |
| 3 | PTN | P21246 | −2.30 (1.4e-21) | −2.14 (1.1e-12) | 3.8% (4/106) q 0.50 | 5 | 0.9830 | 0.0422 | 0.0506 | 0.051 ab | 0.6267 |
| 4 | HSPB6 | O14558 | −2.58 (8.5e-16) | −2.27 (9.8e-16) | 2.8% (3/106) q 0.89 | 7 | 0.9792 | 0.0624 | 0.0091 | 0.062 sm | 0.6133 |
| 5 | BCHE | P06276 | −2.59 (8.8e-22) | −1.75 (9.5e-14) | 6.6% (7/106) q 0.77 | 2 | 0.9779 | 0.7360 | 0.0286 | 0.736 sm | 0.7867 |
| 6 | CFL2 | Q9Y281 | −1.55 (1.0e-11) | −2.20 (1.7e-15) | 3.8% (4/106) q 0.49 | 3 | 0.9717 | 0.5245 | 0.0063 | 0.525 sm | 0.6933 |
| 7 | GREM2 | Q9H772 | −3.55 (5.2e-23) | −1.35 (5.1e-09) | 2.8% (3/106) q 0.93 | 9 | 0.9714 | 0.0214 | 0.0168 | 0.021 sm | 0.4400 |
| 8 | PRPH | P41219 | −1.65 (7.0e-22) | −2.26 (9.5e-16) | 5.7% (6/106) q 0.77 | 12 | 0.9703 | 0.0006 | 0.0008 | 0.001 ab | 0.3467 |
| 9 | TSPAN7 | P41732 | −2.69 (3.9e-23) | −1.00 (4.7e-03) | 4.7% (5/106) q 0.45 | 4 | 0.9665 | 0.7536 | 0.0253 | 0.754 sm | 0.6533 |
| 10 | B3GALT5 | Q9Y2C3 | −2.66 (2.4e-18) | −1.38 (9.5e-09) | 3.8% (4/106) q 1.00 | 10 | 0.9648 | 0.0403 | 0.0987 | 0.099 ab | 0.4267 |
| 11 | CR2 | P20023 | −1.97 (1.7e-13) | −1.41 (2.4e-06) | 9.4% (10/106) q 0.86 | 8 | 0.9646 | 0.1421 | 0.6631 | 0.663 ab | 0.5200 |
| 12 | COLEC12 | Q5KU26 | −1.52 (2.4e-16) | −1.87 (9.3e-14) | 7.5% (8/106) q 0.77 | 11 | 0.9629 | 0.0102 | 0.1310 | 0.131 ab | 0.4000 |
| 13 | MRGPRF | Q96AM1 | −1.74 (6.4e-12) | −1.39 (1.4e-15) | 4.7% (5/106) q 0.76 | 14 | 0.9625 | 0.0096 | 0.0687 | 0.069 ab | 0.2800 |
| 14 | CP | P00450 | −1.68 (7.5e-19) | −1.18 (7.2e-14) | 11.3% (12/106) q 0.39 | 13 | 0.9621 | 0.4024 | 0.2074 | 0.402 sm | 0.3467 |
| 15 | GNAO1 | P09471 | −1.77 (1.3e-19) | −2.23 (9.5e-16) | 3.8% (4/106) q 1.00 | 15 | 0.9610 | 0.1030 | 0.0002 | 0.103 sm | 0.2533 |

All 15: `n_layers_supporting = 2`, concordance "concordant (both down)", `druggable (Y/N) = NA (Lab 3 set no threshold)`, `structure_available` Y only for CA2 and CR2. No CNV data anywhere ("CNV NA").

### Ranking logic
- **Lab 2 score:** equal-weight sum of rank percentiles of |RNA log2FC|, |protein log2FC| and genomic −log10 p, over the 4,851 genes surviving the three-way inner join. No significance/fold-change pre-filter.
- **Lab 3 re-rank:** `final = 0.6 × pct(Module 2 score) + 0.4 × pct(max(p_sm, p_ab))`. Priors = StandardScaler + logistic regression (C = 0.05) on ESM C 300M embeddings vs Open Targets clinical-precedence labels; family-split AUROC 0.751 / AP 0.578 (chance 0.333); random split AUROC 0.833 (reported but rejected as inflated). CA2 removed from the 1,201-protein training panel before fitting. **No druggable/not-druggable threshold.**
- **Deliverable 1 shortlist ranking:** explicitly **lexicographic**, "the same order the evidence was produced" — integration rank first (CA2 1, CR2 11 of 4,851), then the combined omics+tractability score (CA2 1.000, CR2 0.520; Lab 3 ranks 1 and 8 of 15), then pocket confidence. No tie-break needed because CA2 leads on every criterion.
- **"Good target" definition** (replacing the brief's, which no candidate meets): (i) significant and direction-concordant on RNA and protein; (ii) a tractability prior naming a modality; (iii) a confidently predicted structure whose pocket/epitope has its own mean pLDDT ≥ 70. Genomic support dropped from ranking to a shared weakness because *no* candidate has it.
- **Why only 2:** only CA2 and CR2 have all three evidence types; adding TSPAN7/BCHE would mean ranking proteins without the structural check that just rejected CR2's best pocket. Decision recorded in `inventory.md` ("Decisions answered 29 October 2026", item 3).

### The direction problem
- Lab 2 score ranks |log2FC|, so it rewards loss: **43 of the top 50 genes, and all 15 triage candidates, are lower in tumour** than in normal mucosa (verified now from `modules/module2/targets_coad.csv` top 50: 43 down/down, 5 up/up, 2 discordant).
- §2: "my candidates are mostly markers of lost colon identity, which any plan to inhibit them has to deal with. I still took the top 15 forward, because that is what the method produced and re-ranking by hand would have meant picking the answers I wanted."
- §5 CA2: "Its weakness is biological, not structural: CA2 is strongly lost in tumours, so inhibiting it wouldn't undo anything the tumour gained. I think the honest use for CA2 is as a differentiation biomarker, or as a test of whether losing it actually helps the tumour grow." CR2: "Like CA2, CR2 is lower in tumour, so an antibody could only reach the cells that still express it."
- Handoff (DOCX version only): "because my score ignored direction, the **5 genes in the top 50 of `targets_coad.csv` that go up in tumour on both layers (PLAU, TMEM97, SOX9, DHCR7 and REEP6)** should be triaged next, since a protein the tumour gains is a much more natural thing to inhibit." (Verified: top-50 up/up = PLAU rank 22, TMEM97 27, SOX9 40, DHCR7 44, REEP6 49.)
- §6 limitations: "Association, not causation"; "Coverage and direction. Only 2 of 15 candidates had a structural check, and every candidate is lower in tumour."
- Validation plan for CA2: (1) IHC on matched COAD/normal mucosa + qPCR/Western in cell lines, to separate epithelial loss from tissue composition; (2) CRISPR KO / siRNA knockdown in lines retaining CA2 and re-expression in lines that lost it, cross-checked against DepMap — "If CA2 loss is not functionally relevant, the target hypothesis should be dropped in favour of a biomarker use"; (3) DSF thermal shift with acetazolamide as positive control to confirm the assay reads pocket 2.

---

## 11. Disagreements between report text and underlying data/outputs (listed, not reconciled)

**Lab 7**

1. **Working notebook ≠ submitted artifacts.** The current `BIOT6900_Module7_Starter.ipynb` on disk is uncommitted and reverted: cells 3 and 7 fail (no Neo4j), cells 10–18 carry no outputs, and it has **no `CONFIDENCE` dict**, so its `add_edge` cannot produce the `confidence` column that `edges.csv` actually contains, and its cell 12 Cypher does not emit confidence. The README says "the graded notebook, Restart & Run All with outputs". Only the committed copy (`bf72daa`) matches `edges.csv`, `my_knowledge_graph.png` and the report.
2. **Report §1 vs Lab 2/3 ranking: "My strongest RNA-and-protein genes … CA2, DES, CR2, MYH11."** CA2, DES and CR2 are Lab 2 ranks 1, 2 and 11 and are on the 15-gene triage list, but **MYH11 is Lab 2 rank 24 and is not on that list**; the genes at ranks 3–10 (PTN, HSPB6, BCHE, CFL2, GREM2, PRPH, TSPAN7, B3GALT5) were skipped without comment.
3. **Report §2 Table 2 omits two genes present in the PageRank output.** JAK2 (0.0158) and IL23R (0.0086) appear in the cell-14 output but not in the report's table, which is presented as the biomarker ranking.
4. **Report's path string drops the confidence tokens.** Cell 12 printed `TARGETS [ChEMBL, high] … HAS_PHENOTYPE [HPO, low]`; §3 renders them as `[ChEMBL]`, `[GWAS Catalog, p = 1 × 10⁻²¹]`, `[HPO]`, `[HPO]` — a p-value substituted for a confidence label, with no edge in `edges.csv` carrying a pvalue on the TARGETS/HAS_PHENOTYPE hops.
5. **Report §2 says community detection put "DES, MYH11 and CR2 in one cluster … and two symptoms"** — correct, but the same community also contains the COAD disease node itself (stated) while CA2 is in a *separate* community (9). Any reading of community 1 as "the tissue-composition cluster" excludes CA2 even though CA2 is the clearest "lost in tumour" gene in the dataset and is described the same way in D1.
6. **"PageRank fused with your gene evidence" (notebook cell 14 header, report Table 2 caption)** — nothing is actually fused; `pagerank`, `gwas_p` and `druggability` are printed side by side. The report itself says so in §2 ("My table shows PageRank and evidence side by side instead of merging them into one number"), so the caption and the body disagree.
7. **p-value float noise.** `nodes.csv`/`edges.csv` store MYC as `1.9999999999999998e-122` and SMAD7 as `3.0000000000000003e-30`, while `gwas_lookup.csv`/`targets_coad.csv` and the report give 2e-122 and 3e-30. Cosmetic, but the graph artifacts do not literally contain the reported values.
8. **Report §4 "Phenotype edges (8). … all five COAD phenotypes are very general"** — 8 is the total including the 3 UC edges; only 5 are COAD's. Not wrong, but the count and the clause refer to different sets.
9. **rs7229639 absent.** The Module 6 PRS locus at the SMAD7 region (chr18:46,452,637) was dropped in Module 6 for a missing effect allele and never reached the graph; the report's claim that the variant layer uses "four COAD GWAS lead SNPs from Module 6" is consistent with the files, but a reader expecting the full Module 6 SMAD7 locus set will not find rs7229639.
10. **Repurposing drug status.** Ustekinumab and Tofacitinib carry `status = "approved (UC)"` in `nodes.csv`, yet the 2-hop COAD query output ("approved (CRC)" / "approved (other indication)") is what §3 narrates; the repurposing candidates are never described with their stored status.

**Deliverable 1**

11. **CR2 vs 1LY2 RMSD:** `modules/module4/structure_assessment_coad.csv` records **38.45 Å**; D1 §4 uses **2.20 Å over 80 residues** for Sushi domain 1 and calls the saved value an alignment artefact. `inventory.md` gap 7 notes the earlier side-check figure was **1.00 Å over residues 20–149** and "isn't saved to any file" — so three different numbers circulate (38.45 saved, 2.20 in D1 text per `provenance.md` NB4 cell 12, 1.00 in `inventory.md` / the Lab 4 report).
12. **CR2 Sushi-repeat count:** `deliverable1.md` §4 says PAE splits the extracellular region into "**eight confident domains** matching CR2's Sushi repeats"; the final DOCX says "eight confident blocks that follow CR2's **15** Sushi repeats". Eight blocks vs 15 repeats is not the same claim.
13. **CR2 antibody-prior comparison:** `deliverable1.md` §3 says "the next-highest antibody prior is COLEC12 at 0.131"; the DOCX says "the next largest antibody-over-small-molecule **gap** is COLEC12 at +0.121". Different statistics, same sentence slot.
14. **Figure 1B region count:** the caption says the 15 triage candidates "all sit in the RNA ∩ protein **only** region (3,099 genes)", while §2 text reports 3,115 for RNA ∩ protein; `provenance.md` resolves them as different quantities (3,115 = RNA∩P, 3,099 = RNA∩P∧¬mutation), but the two numbers appear one paragraph apart with no distinction drawn in the prose. `inventory.md` also records RNA∩protein = 3,115 and triple = 16, i.e. 3,099 = 3,115 − 16.
15. **Genomic FDR denominator:** `inventory.md` gap 10 — the Module 2 report says "69 genes significant at FDR < 0.05" on the genomic layer (all 14,706 genes), while D1 reports **26** (within the 4,851 joined genes). Both documents stand.
16. **"Druggable (Y/N)" column** in `evidence_table.csv` is `NA (Lab 3 set no threshold)` for all 15, yet `inventory.md` Decision 2 option (b) lists what a ≥ 0.5 cut-off would give (CA2 0.770, TSPAN7 0.754, BCHE 0.736, CR2 0.663, CFL2 0.525 = Y). The published table reports no threshold; the inventory records the rejected alternative.
17. **Brief vs delivery:** the brief asked for a 3–5 target shortlist and for "the triple-overlap region contains the strongest candidates"; D1 delivers 2 targets and states the triple-overlap framing is **false for this dataset** (the 16 triple-overlap genes contain none of the 15 candidates). Deliberate and documented, but a disagreement with the assignment text.
18. **Two reference entries are flagged `[VERIFY]`** in `deliverable1.md` (EvolutionaryScale 2024 ESM Cambrian; Ochoa et al. 2023 Open Targets) — unresolved at the time the markdown was written.

**Cross-module (Lab 7 vs D1)**

19. **Same genes, opposite framings of CA2/CR2.** D1 ranks CA2 #1 and CR2 #2 as the lead target package; Lab 7's graph puts CA2 mid-table (PageRank 0.0158) and CR2 at the floor (0.0086), and Lab 7 §2 reads the DES/MYH11/CR2 cluster as a tissue-composition artefact — "a warning about my own Module 2 ranking", i.e. a warning about the list D1 is built on. Neither document cites the other.
20. **Druggability values are from different fits.** `prepare_kg_data.py` *refits* the Module 3 prior with every graph gene held out, so Lab 7's `druggability` for CA2 is **0.769118** vs Module 3's published `p_sm` **0.769844** (and CR2 0.64031 vs 0.663090). The script prints a "refit prior vs Module 3" comparison, but the Lab 7 report quotes the values as if they were the Module 3 prior ("my Module 3 tractability prior, max(p_sm, p_ab)"), and D1 Table 1 quotes the Module 3 numbers (0.770 / 0.663). Same label, two sets of numbers.
21. **DES druggability:** `nodes.csv` 0.0015295037 (= refit `p_ab`); Module 3 `p_sm` 0.000992 / `p_ab` 0.001368; D1 reports 0.001. Lab 7 Table 2 says the bottom group spans "0.00–0.64", consistent, but no single source gives DES the same value twice.

---

## MISSING items

- `[MISSING: rs7229639 anywhere in Module 7]` — see §9. Present only in Module 6 (dropped from the PRS for a missing effect allele).
- `[MISSING: CA2 in the DES/MYH11/CR2 community]` — CA2 is in community 9, not community 1.
- `[MISSING: working-copy notebook outputs for cells 10–18]` — no saved outputs in the on-disk notebook; all query/PageRank/community/figure outputs were read from the committed version (`bf72daa`).
- `[MISSING: Part I data under /Users/drashtibhanushali/biot6900/shared/module2_omics/ beyond the per-layer TSVs]` — the folder holds `PROVENANCE.md`, `coad_transcriptomics.tsv`, `coad_proteomics.tsv`, `coad_genomics.tsv`, `coad_kras_variant_spectrum.tsv`, `coad_cptac_matched_rna_protein_corr.tsv`, `coad_independent_proteomics_fudan.tsv`, three CPTAC-BRCA files and `raw/`; no Lab 7-specific data.
- `[MISSING: a saved file holding the CR2 1.00 Å side-check RMSD]` — per `inventory.md` gap 7, that value exists only in prose.
- `[MISSING: nodes.csv gwas_p for the UC genes' per-SNP p-values on the OCCURS_IN edges]` — UC lead SNPs (rs11209026, rs10758669, rs6871626) are in `gwas_lookup.csv` but are **not** variant nodes in the graph; only the four COAD lead SNPs are.
