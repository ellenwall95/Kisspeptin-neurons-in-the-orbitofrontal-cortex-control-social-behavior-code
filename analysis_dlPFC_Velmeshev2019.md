# Analysis of KISS1 and KISS1R Expression in the Human Dorsolateral Prefrontal Cortex (dlPFC)

**Dataset:** Velmeshev et al. (2019) *Science*, 364, 685
Brodmann areas 9/46 · 23 donors (10 Control, 13 ASD) · PFC samples only (ACC/BA24 excluded)

**Figures produced:**
- Main Figure 6i — Cell-type enrichment (all donors pooled)
- Main Figure 6j — UMAP with KISS1/KISS1R expression overlay
- Supplementary Fig. 6e — KISS1 per-donor % by diagnosis (Control vs. ASD)
- Supplementary Fig. 6f — KISS1R per-donor % by diagnosis (Control vs. ASD)

**Software:** Python 3.12 · pandas 2.1 · scipy 1.12 · statsmodels 0.14

---

## Section 1 — Configuration

| Parameter | Value |
|---|---|
| Expression matrix | `velmeshev2019_expression_matrix.csv` |
| Cell metadata | `velmeshev2019_cell_metadata.csv` (Data S2) |
| UMAP coordinates | `velmeshev2019_umap_coords.csv` (if separate) |
| Minimum nuclei per donor × cell-type | 5 |
| Output directory | `results_dlPFC/` |

**Neuronal cell types (11):**
L5/6, L5/6-CC, L2/3, L4, IN-PV, IN-SST, IN-SV2C, IN-VIP, Neu-NRGN-I, Neu-NRGN-II, Neu-mat

---

## Section 2 — Data Loading

Load expression matrix and cell-level metadata. Two supported formats:

- **Option A (h5ad):** Use `anndata.read_h5ad()` and extract as needed.
- **Option B (CSV/tabular):** Read metadata and expression matrix separately; look up KISS1 and KISS1R rows by gene symbol.

Build one row per nucleus with columns: `donor`, `cell_type`, `diagnosis`, `sex`, `region`, `KISS1`, `KISS1R`.

**Filter:** Retain only `region == "PFC"` (excludes ACC / BA24). Then restrict to the 11 neuronal cell types.

---

## Section 3 — Cell-Type Enrichment (Fig. 6i)

All 23 PFC donors pooled regardless of diagnosis.

> **Note on sex stratification:** Sex-stratified analysis was not performed for this dataset. Only 3 female neurotypical controls contributed PFC samples, with zero KISS1-expressing nuclei and only 2 KISS1R-expressing nuclei detected across all female control cells.

**Method:** One-vs-rest Fisher's exact test (two-sided) for each cell type × gene.

**Multiple testing:** BH correction across 11 cell types, applied separately per gene.

**Output columns:** `gene`, `cell_type`, `n_positive`, `n_total`, `pct_expressing`, `odds_ratio`, `pval`, `qval`, `significant`

**Output file:** `results_dlPFC/Fig6i_celltype_enrichment.csv`

### Contingency table per cell type

|  | Gene+ | Gene− |
|---|---|---|
| **In cell type** | pos_in | neg_in |
| **Outside cell type** | pos_out | neg_out |

---

## Section 4 — UMAP Visualization (Fig. 6j)

Requires `UMAP_1` and `UMAP_2` columns in the neuronal dataframe. If absent, this section is skipped with a warning.

**Three panels:**

| Panel | Content | Color |
|---|---|---|
| 1 | Cell-type reference (tab20 palette) | Per cell type |
| 2 | KISS1-expressing nuclei highlighted | steelblue |
| 3 | KISS1R-expressing nuclei highlighted | mediumseagreen |

Non-expressing nuclei rendered in light grey (size 0.5, alpha 0.1). Expressing nuclei rendered at size 3, alpha 0.7.

**Output file:** `results_dlPFC/Fig6j_UMAP_KISS1_KISS1R.pdf` (300 dpi)

---

## Section 5 — Per-Donor Percentages

For every donor × cell type × diagnosis combination with ≥ 5 nuclei, compute:

$$\text{pct} = \frac{\text{nuclei with gene} > 0}{\text{total nuclei}} \times 100$$

Donors with PFC nuclei but zero expression are assigned 0%. Donor × cell-type combinations below the 5-nucleus threshold are excluded.

**Output columns:** `donor`, `cell_type`, `diagnosis`, `gene`, `n_total`, `n_positive`, `pct`

**Output file:** `results_dlPFC/per_donor_pct_all.csv`

---

## Section 6 — Control vs. ASD Comparison (Suppl. Fig. 6e–f)

Per-donor percentages (from Section 5) compared between neurotypical Controls and ASD cases using a **two-sided Mann–Whitney U test** for each gene × cell type.

**Multiple testing:** BH correction across 11 cell types, applied separately per gene.
NaN p-values (insufficient donors in one group) are excluded from correction and remain NaN in output.

**Output columns:** `gene`, `cell_type`, `ctrl_n`, `asd_n`, `ctrl_median_pct`, `asd_median_pct`, `U_stat`, `pval`, `qval`, `significant`

**Output file:** `results_dlPFC/SuppFig6ef_control_vs_ASD.csv`

---

## Execution Order

```
load_data()
  └─ filter_neuronal()
        ├─ celltype_enrichment()   →  Fig6i_celltype_enrichment.csv
        ├─ plot_umap()             →  Fig6j_UMAP_KISS1_KISS1R.pdf
        ├─ per_donor_pct()         →  per_donor_pct_all.csv
        └─ control_vs_asd()        →  SuppFig6ef_control_vs_ASD.csv
```
