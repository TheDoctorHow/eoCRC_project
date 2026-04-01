# eoCRC Gut Microbiome Meta-Analysis

This pipeline integrates nine international public cohorts accessed through the curatedMetagenomicData v3 Bioconductor package, applies MaAsLin2 differential abundance testing with proper centered log-ratio (CLR) pre-computation, and uses Leave-One-Study-Out cross-validation (LOSO-CV) to obtain honest cross-cohort generalization estimates. Random forest classification, XGBoost, and ElasticNet models are evaluated; SHAP values provide feature-level interpretability. The primary finding is a reproducible eoCRC microbiome signal dominated by oral pathobionts and depleted butyrate producers, with species-level AUROC of 0.704 (95% CI 0.602–0.789, perm-p = 0.002) and low between-cohort heterogeneity (I² = 4.5%).

---

## Dataset

| Group | Definition | N |
|---|---|---|
| eoCRC | CRC, age < 50 | 78 |
| loCRC | CRC, age >= 50 | 562 |
| Young controls | Healthy, age < 50 | 121 |
| Older controls | Healthy, age >= 50 | 521 |
| **Total** | | **1,282** |

**Cohorts (9):** ZellerG_2014, FengQ_2015, YuJ_2015, WirbelJ_2018, ThomasAM_2019_c, VogtmannE_2016, YachidaS_2019, GuptaA_2019, HanniganGD_2017

**Feature matrices (post-prevalence filter, >= 10%):** 220 species, 419 unstratified pathways

**Source:** curatedMetagenomicData v3 (Bioconductor) — MetaPhlAn3 taxonomic + HUMAnN3 functional profiles, uniformly processed.

---

## Key Results

| Analysis | Result |
|---|---|
| Primary AUROC (species RF, LOSO-CV) | **0.704** (95% CI 0.602–0.789) |
| Between-cohort heterogeneity | I² = 4.5% |
| Permutation p-value | 0.002 |
| Significant species (MaAsLin2 v3, q < 0.25, eoCRC vs young controls) | **47** (26 enriched, 21 depleted) |
| Doubly validated SHAP features (SHAP top-20 AND MaAsLin2 q < 0.25) | **13 / 20** |
| Top 3 SHAP species | Parvimonas micra, Peptostreptococcus stomatis, Gemella morbillorum |

**Key biological signal:** Enrichment of oral pathobionts (Parvimonas micra, Gemella morbillorum, Peptostreptococcus stomatis, Dialister pneumosintes, Anaerotruncus colihominis) and depletion of butyrate producers (Faecalibacterium prausnitzii, Lactobacillus rogosae, Turicimonas muris) in eoCRC relative to young healthy controls.

---

## Repository Structure

```
scripts/       R scripts numbered step1-step5 (main pipeline) + figure regeneration scripts
figures/       Output figures (PNG + PDF), 300 dpi
results/       CSV outputs: LOSO-CV summary, SHAP values, MaAsLin2 v3 summaries
README.md      This file
summary.md     Detailed results narrative
.gitignore     Excludes large .rds data objects, old v2 results
```

---

## Software Requirements

| Software | Version |
|---|---|
| R | 4.5.3 |
| Bioconductor | 3.22 |
| curatedMetagenomicData | 3.18.0 |
| SummarizedExperiment | 1.40.0 |
| Maaslin2 | 1.24.1 |
| compositions | 2.0.9 |
| randomForest | 4.7.1.2 |
| xgboost | 3.2.1.1 |
| glmnet | 4.1.10 |
| pROC | 1.19.0.1 |
| meta | 8.2.1 |
| shapviz | 0.10.3 |
| ggplot2 | 4.0.2 |
| ggrepel | 0.9.8 |
| pheatmap | 1.0.13 |
| patchwork | 1.3.2 |
| dplyr | 1.2.0 |
| tidyr | 1.3.2 |

Install all packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Maaslin2", "SummarizedExperiment", "curatedMetagenomicData"))
install.packages(c("compositions", "randomForest", "xgboost", "glmnet",
                   "pROC", "meta", "shapviz", "ggplot2", "ggrepel",
                   "pheatmap", "patchwork", "dplyr", "tidyr"))
```

---

## Reproduction Instructions

1. **Clone the repository**

   ```bash
   git clone https://github.com/TheDoctorHow/eoCRC_pipeline.git
   cd eoCRC_pipeline
   ```

2. **Download raw data** (requires Bioconductor + internet access; not stored in this repo due to file size)

   ```r
   # Run from the project root directory
   Rscript scripts/step1_extract_filter.R
   ```

   This fetches data via `curatedMetagenomicData`, extracts species and pathway matrices, applies prevalence filtering, and saves all `.rds` data objects to the project root.

3. **Run MaAsLin2 differential abundance** (v3 — proper CLR pre-computation)

   ```r
   Rscript scripts/step2c_maaslin2_v3.R
   ```

   Outputs to `maaslin2_results_v3/`. Step2 and step2b are retained for methodological reference but step2c is canonical.

4. **Run LOSO-CV classification**

   ```r
   Rscript scripts/step3_loso_cv.R       # primary: eoCRC vs young controls
   Rscript scripts/step3b_secondary_perm.R  # secondary: eoCRC vs loCRC
   ```

5. **Compute SHAP feature importance**

   ```r
   Rscript scripts/step4_shap.R
   ```

6. **Regenerate all figures**

   ```r
   Rscript scripts/regen_fig1.R
   Rscript scripts/regen_fig2.R
   Rscript scripts/regen_fig3.R
   Rscript scripts/regen_fig4.R
   Rscript scripts/step5_figures.R   # figures 5-7
   ```

   All figures are saved to `figures/`.

---

## Methods Note — CLR Pre-computation

MaAsLin2 v1.24.1's internal CLR normalization uses `chemometrics::clr(x + 1)`. For relative abundance proportions (x << 1), `log(1 + x) ≈ x`, which reduces to arithmetic mean-centering rather than true log-ratio CLR. All results here use externally pre-computed CLR via `compositions::clr(mat + 1e-6)` passed to MaAsLin2 with `normalization = "NONE"` (see `scripts/step2c_maaslin2_v3.R`). See `results/maaslin2_v3_summaries/` for the canonical association results.

