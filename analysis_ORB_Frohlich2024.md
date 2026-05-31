# Analysis of KISS1 and KISS1R Expression in the Human Orbitofrontal Cortex (OFC)

**Dataset:** Fröhlich et al. (2024) *Nature Neuroscience*
Brodmann area 11 · 87 post-mortem donors · GEO accession: GSE254569

**Figures produced:**
- Main Figure 6h — Sex-stratified cell-type enrichment (controls only)
- Supplementary Fig. 6a — KISS1 per-donor % in excitatory neurons by diagnosis
- Supplementary Fig. 6b — KISS1R per-donor % in excitatory neurons by diagnosis
- Supplementary Fig. 6c — KISS1 per-donor % in inhibitory neurons by diagnosis
- Supplementary Fig. 6d — KISS1R per-donor % in inhibitory neurons by diagnosis

**Software:** Python 3.12 · pandas 2.1 · scipy 1.12 · statsmodels 0.14

---

## Section 1 — Configuration

| Parameter | Value |
|---|---|
| H5AD file | `GSE254569_OFC_snRNAseq.h5ad` |
| KISS1 Ensembl ID | `ENSG00000170498` |
| KISS1R Ensembl ID | `ENSG00000116014` |
| Minimum nuclei per donor × cell-type | 5 |
| Output directory | `results_OFC/` |

**Excitatory cell types (8):**
Exc_L2–L3, Exc_L3–L5, Exc_L4–L6_1, Exc_L4–L6_2, Exc_L4–L6_3, Exc_L5–L6_1, Exc_L5–L6_2, Exc_L5–L6_HTR2C

**Inhibitory cell types (7):**
In_PVALB_Ba, In_PVALB_Ch, In_SST, In_VIP, In_RELN, In_LAMP5_1, In_LAMP5_2

**Diagnostic groups:** Control, Schizophrenia, Schizoaffective, Bipolar, MDD

---

## Section 2 — Data Loading

Load the Fröhlich et al. h5ad file and extract KISS1/KISS1R raw counts alongside cell-level metadata.

**Gene lookup order:**
1. Match by Ensembl ID in `adata.var_names`
2. Match by gene symbol in `adata.var_names`
3. Match by gene symbol in `adata.var["gene_name"]` column

Handle sparse matrices (`toarray()` if needed). Build one row per nucleus with columns: `donor`, `cell_type`, `diagnosis`, `sex`, `KISS1`, `KISS1R`.

Filter to the 15 neuronal cell types (excitatory + inhibitory) before any analysis.

---

## Section 3 — Cell-Type Enrichment, Sex-Stratified (Fig. 6h)

Restricted to **neurotypical control donors** only.

**Method:** One-vs-rest Fisher's exact test (two-sided) for each cell type × gene × sex.

**Multiple testing:** Benjamini–Hochberg (BH) correction across 15 cell types, applied separately per gene per sex.

**Output columns:** `gene`, `sex`, `cell_type`, `n_positive`, `n_total`, `pct_expressing`, `odds_ratio`, `pval`, `qval`, `significant`

**Output file:** `results_OFC/Fig6h_celltype_enrichment_sex_stratified.csv`

### Contingency table per cell type

|  | Gene+ | Gene− |
|---|---|---|
| **In cell type** | pos_in | neg_in |
| **Outside cell type** | pos_out | neg_out |

---

## Section 4 — Per-Donor Percentages

For every donor × cell type × diagnosis × sex combination with ≥ 5 nuclei, compute:

$$\text{pct} = \frac{\text{nuclei with gene} > 0}{\text{total nuclei}} \times 100$$

Donor × cell-type combinations below the 5-nucleus threshold are excluded entirely (not assigned 0%).

**Output columns:** `donor`, `cell_type`, `diagnosis`, `sex`, `gene`, `n_total`, `n_positive`, `pct`

**Output file:** `results_OFC/per_donor_pct_all.csv`

---

## Section 5 — Diagnostic Group Comparisons (Suppl. Fig. 6a–d)

Each non-Control diagnostic group (Schizophrenia, Schizoaffective, Bipolar, MDD) is compared to Controls in a **two-stage test** for each gene × cell type combination.

### Stage 1 — Pooled cell-level Fisher's exact test

All nuclei from the diagnostic group vs. all Control nuclei, two-sided.

### Stage 2 — Donor-level Mann–Whitney U test

Per-donor percentages (from Section 4) compared between groups, two-sided.

**Multiple testing:** BH correction applied across 15 cell types separately for Fisher p-values and MWU p-values.
NaN p-values (groups with < 1 donor) are excluded from correction and remain NaN in output.

**Output columns:** `gene`, `diagnosis`, `cell_type`, `fisher_pval`, `fisher_qval`, `mwu_pval`, `mwu_qval`, `ctrl_n_donors`, `dx_n_donors`, `ctrl_median_pct`, `dx_median_pct`

**Output file:** `results_OFC/SuppFig6a-d_diagnostic_comparisons.csv`

---

## Execution Order

```
load_data()
  └─ filter_neuronal()
        ├─ celltype_enrichment_by_sex()    →  Fig6h_celltype_enrichment_sex_stratified.csv
        ├─ per_donor_pct()                 →  per_donor_pct_all.csv
        └─ diagnostic_group_comparisons()  →  SuppFig6a-d_diagnostic_comparisons.csv
```
