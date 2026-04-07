# Reproducing the eoCRC Microbiome Meta-Analysis

Step-by-step instructions to reproduce the full pipeline on a fresh machine.

---

## Prerequisites

### R version
R ≥ 4.2.0 (tested on R 4.4.x)

### Bioconductor packages
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "curatedMetagenomicData",
  "TreeSummarizedExperiment",
  "SummarizedExperiment",
  "Maaslin2"
))
```

### CRAN packages
```r
install.packages(c(
  "compositions", "randomForest", "xgboost", "glmnet",
  "pROC", "meta", "fastshap", "shapviz",
  "ggplot2", "ggrepel", "patchwork", "pheatmap",
  "dplyr", "tidyr", "parallel", "RColorBrewer"
))
```

> **Note:** `xgboost` ≥ 1.7 uses the v3 API (`xgb.DMatrix` / `xgb.train`). The scripts are written for this API. Older versions will fail.

---

## Cloning the repository

```bash
git clone https://github.com/TheDoctorHow/eoCRC_pipeline.git
cd eoCRC_pipeline
```

All scripts must be run **from the repo root** (`eoCRC_pipeline/`), not from inside `scripts/`. Each script auto-detects its working directory:

```r
if (basename(getwd()) == "scripts") setwd("..")
```

---

## Running the pipeline

### Step 1 — Data download and prevalence filtering

```bash
Rscript scripts/step1_extract_filter.R
```

**What it does:** Downloads relative abundance and pathway abundance profiles for 9 CRC cohorts from curatedMetagenomicData v3 (Bioconductor), applies a ≥10% prevalence filter, and saves filtered matrices.

**Expected outputs:**

| File | Description | Approx. size |
|------|-------------|-------------|
| `tse_relabund_full.rds` | Full species TSE (908 sp × 1282 samples) | ~1.5 MB |
| `tse_pathways_full.rds` | Full pathway TSE (37,945 pw × 1282 samples) | ~36 MB |
| `all_samples_metadata.rds` | Sample metadata (1282 rows × 6 cols) | <1 MB |
| `species_filtered.rds` | Filtered species matrix (1282 × 220) | ~640 KB |
| `path_filtered.rds` | Filtered pathway matrix (1282 × 419) | ~3 MB |
| `species_filtered.csv` | Same as above, CSV format | ~1.2 MB |
| `path_filtered.csv` | Same as above, CSV format | ~5 MB |
| `species_prevalence.rds` | Per-species prevalence vector | <1 MB |
| `path_prevalence.rds` | Per-pathway prevalence vector | <1 MB |

**Expected warnings during download:** Messages like `dropping rows without rowTree matches: k__Bacteria|...` are printed by curatedMetagenomicData internally for a handful of taxa per cohort that lack phylogenetic tree entries. These are harmless and expected — those taxa are still included in the abundance matrix.

**Re-running:** If RDS files already exist on disk, step1 loads from cache instead of re-downloading. Delete the three `tse_*.rds` and `all_samples_metadata.rds` files to force a fresh download.

---

### Step 2c — MaAsLin2 differential abundance (canonical version)

```bash
Rscript scripts/step2c_maaslin2_v3.R
```

**What it does:** Runs MaAsLin2 with pre-computed CLR (using `compositions::clr`, pseudocount=1e-6) for three comparisons: eoCRC vs Young Controls, loCRC vs Older Controls, eoCRC vs loCRC. Both species and pathway feature sets.

**Requires:** Output from Step 1.

**Expected outputs:**

| File/Dir | Description |
|----------|-------------|
| `maaslin2_results_v3/` | Per-comparison MaAsLin2 output directories |
| `maaslin2_results_v3/summaries/all_associations.csv` | All associations, all comparisons |
| `maaslin2_results_v3/summaries/significant_q0.25.csv` | Significant hits (q < 0.25) |

**Key result:** 47 significant species in eoCRC vs Young Controls (q < 0.25).

> Steps 2 and 2b (`step2_maaslin2.R`, `step2b_maaslin2_corrected.R`) are earlier versions kept for reference. Step 2c is the canonical run used for all figures.

---

### Step 3 — LOSO-CV machine learning classification

```bash
Rscript scripts/step3_loso_cv.R
```

**What it does:** Leave-One-Study-Out cross-validation with three models (Random Forest, XGBoost, ElasticNet), three feature sets (species, pathways, combined), permutation test (1000 shuffles), and learning curves. Primary: eoCRC vs Young Controls. Secondary: eoCRC vs loCRC.

**Requires:** Output from Step 1.

**Expected outputs:**

| File | Description |
|------|-------------|
| `loso_primary_folds.rds` | Per-fold AUROC results, primary analysis |
| `loso_secondary_folds.rds` | Per-fold AUROC results, secondary analysis |
| `loso_primary_permutations.rds` | 1000 permutation null AUROCs, primary |
| `loso_secondary_permutations.rds` | 1000 permutation null AUROCs, secondary |
| `loso_primary_learning_curves.rds` | Learning curve data |
| `loso_pooled_auroc_summary.rds` | D-L pooled AUROC table |
| `loso_pooled_auroc_summary.csv` | Same, CSV format |

**Key result:** Species RF pooled AUROC = 0.704 (95% CI 0.602–0.789), I² = 4.5%, perm-p = 0.002.

> **Note:** Step 3b (`step3b_secondary_perm.R`) and Step 3c (`step3c_xgb_patch.R`) are patch scripts that were used during development to fix specific issues. They are not needed in a fresh run of step3.

---

### Step 4 — SHAP feature interpretation

```bash
Rscript scripts/step4_shap.R
```

**What it does:** Trains a final Random Forest on all 199 eoCRC + Young Control samples, computes approximate SHAP values (fastshap, nsim=100), cross-references top-20 SHAP features with MaAsLin2 results.

**Requires:** Output from Step 1. (MaAsLin2 results are read from `maaslin2_results_v2/` for the cross-reference step — run Step 2b if you need this cross-reference, otherwise the main SHAP computation still runs.)

**Expected outputs:**

| File | Description |
|------|-------------|
| `shap_values.rds` / `.csv` | Full SHAP matrix (199 samples × 50 features) |
| `shap_top20_summary.rds` / `.csv` | Top-20 feature importance table |
| `figures/fig_shap_beeswarm.png` | Beeswarm importance plot |
| `figures/fig_shap_bar.png` | Bar importance plot |
| `figures/fig_shap_dependence.png` | Dependence plots (top 5) |
| `figures/fig_shap_combined.png` | Combined beeswarm + bar |

**Key result:** 13 doubly-validated features (in both SHAP top-20 and MaAsLin2 q < 0.25).

---

### Step 5 — All presentation figures

```bash
Rscript scripts/step5_figures.R
```

**What it does:** Generates all 7 presentation figures plus a `summary.md` report.

**Requires:** Outputs from Steps 1, 3, and 4.

**Expected outputs in `figures/`:** All PNG (300 dpi) + PDF pairs for figures 1–7.

---

### Regenerating individual figures

Each figure can be regenerated independently without re-running the full pipeline:

```bash
Rscript scripts/regen_fig1.R   # Study design / LOSO schematic
Rscript scripts/regen_fig2.R   # AUROC forest plot (hardcoded data, no dependencies)
Rscript scripts/regen_fig3.R   # Volcano plot (requires step2c output)
Rscript scripts/regen_fig4.R   # SHAP beeswarm (requires step4 output)
Rscript scripts/regen_fig5.R   # Species boxplots (requires step1 + step2c output)
Rscript scripts/regen_fig6.R   # Pathway heatmap (requires step1 + step2b + step2c output)
Rscript scripts/regen_fig7.R   # Biological narrative (no data dependencies)
Rscript scripts/regen_shap_supp.R  # SHAP supplementary figures (requires step4 output)
```

---

## Recommended run order

```
step1  →  step2c  →  step3  →  step4  →  step5
                  ↘  step2b (needed for fig6 and gen_summary only)
```

Minimum to reproduce main results and figures: **step1 → step2c → step3 → step4 → step5**

---

## Possible errors and fixes

| Error message | Cause | Fix |
|---|---|---|
| `Error: Run scripts/step1_extract_filter.R first` | Missing filtered data files | Run step1 |
| `Error: Run scripts/step3_loso_cv.R first` | Missing LOSO fold RDS files | Run step3 |
| `Error: Run scripts/step4_shap.R first` | Missing `shap_extended_results.rds` | Run step4 |
| `Error in gzfile(file, "rb") : cannot open the connection` | Old script trying to load non-existent RDS | Update to latest scripts from this repo |
| `Error: package 'curatedMetagenomicData' not found` | Bioconductor package not installed | Run `BiocManager::install("curatedMetagenomicData")` |
| XGBoost column name error | xgboost < 1.7 (old API) | Upgrade: `install.packages("xgboost")` |
| `cairo_pdf` not found | Cairo graphics not available | Install: `install.packages("Cairo")` or use `pdf()` fallback |
| Step1 very slow (>30 min) | curatedMetagenomicData downloads each cohort separately | Normal for first run; subsequent runs use cache |


Step 3 is the bottleneck. It uses `parallel::mclapply` with up to 4 cores automatically — running on a machine with ≥4 cores cuts permutation time roughly by 4×.

---
