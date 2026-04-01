#!/usr/bin/env Rscript
# Step 2: MaAsLin2 differential abundance — three group comparisons × two feature sets
# Groups:  eoCRC (CRC, age<50) | loCRC (CRC, age≥60) | Control (all controls)
# Comparisons: eoCRC vs Control | loCRC vs Control | eoCRC vs loCRC
# Features: species_filtered (÷100 → CLR) and path_filtered (CLR)

suppressPackageStartupMessages({
  library(Maaslin2)
  library(dplyr)
  library(tidyr)
})

setwd("/home/yugiu/eoCRC_analysis")

# ── 1. Load data ──────────────────────────────────────────────────────────────
cat("Loading data...\n")
meta              <- readRDS("all_samples_metadata.rds")
species_filtered  <- readRDS("species_filtered.rds")
path_filtered     <- readRDS("path_filtered.rds")

# Species matrix: values are % relative abundance → divide by 100
species_prop <- species_filtered / 100

cat(sprintf("  metadata:  %d x %d\n", nrow(meta), ncol(meta)))
cat(sprintf("  species:   %d x %d  (after /100: range %.2e – %.2e)\n",
            nrow(species_prop), ncol(species_prop),
            min(species_prop), max(species_prop)))
cat(sprintf("  pathways:  %d x %d  (range %.2e – %.2e)\n",
            nrow(path_filtered), ncol(path_filtered),
            min(path_filtered), max(path_filtered)))

# ── 2. Define groups ──────────────────────────────────────────────────────────
meta <- meta %>%
  mutate(
    maaslin_group = case_when(
      group == "CRC"     & age < 50  ~ "eoCRC",
      group == "CRC"     & age >= 60 ~ "loCRC",
      group == "control"             ~ "Control",
      TRUE                           ~ NA_character_
    )
  )

group_counts <- table(meta$maaslin_group, useNA = "ifany")
cat("\nGroup definitions (age gap 50–59 excluded from CRC):\n")
print(group_counts)

# Confirm row alignment between matrices and metadata
stopifnot(all(rownames(species_prop) == meta$sample_id))
stopifnot(all(rownames(path_filtered) == meta$sample_id))
cat("\nRow alignment check passed.\n\n")

# ── 3. Helper: run one MaAsLin2 comparison ───────────────────────────────────
run_maaslin <- function(feat_mat, meta_sub, case_label, ref_label,
                        output_dir, feat_type, covariates) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Build per-run metadata
  df_meta <- meta_sub %>%
    select(sample_id, study_name, maaslin_group, BMI, age) %>%
    mutate(is_case = as.numeric(maaslin_group == case_label)) %>%
    as.data.frame()
  rownames(df_meta) <- df_meta$sample_id
  df_meta$sample_id <- NULL

  # Subset feature matrix
  feat_sub <- feat_mat[rownames(df_meta), , drop = FALSE]

  # Build fixed/random effects
  fixed_effects  <- c("is_case", covariates)
  random_effects <- "study_name"

  cat(sprintf("  Running MaAsLin2: %s vs %s (%s) — %d samples, %d features\n",
              case_label, ref_label, feat_type,
              nrow(feat_sub), ncol(feat_sub)))
  cat(sprintf("  Fixed effects: %s | Random: %s\n",
              paste(fixed_effects, collapse=", "), random_effects))

  fit <- Maaslin2(
    input_data      = feat_sub,
    input_metadata  = df_meta,
    output          = output_dir,
    fixed_effects   = fixed_effects,
    random_effects  = random_effects,
    normalization   = "CLR",
    transform       = "NONE",
    analysis_method = "LM",
    max_significance = 0.25,   # q-value threshold for reporting
    correction      = "BH",
    min_prevalence  = 0,       # already filtered upstream
    min_abundance   = 0,
    cores           = 4,
    plot_heatmap    = FALSE,
    plot_scatter    = FALSE
  )
  invisible(fit)
}

# ── 4. Define three comparisons ───────────────────────────────────────────────
comparisons <- list(
  list(case = "eoCRC",   ref = "Control", tag = "eoCRC_vs_Control"),
  list(case = "loCRC",   ref = "Control", tag = "loCRC_vs_Control"),
  list(case = "eoCRC",   ref = "loCRC",   tag = "eoCRC_vs_loCRC")
)

# Covariates: BMI has 10 NAs — MaAsLin2 drops those rows automatically
# age is included only in eoCRC vs loCRC to guard against residual age confounding
covariate_map <- list(
  eoCRC_vs_Control = "BMI",
  loCRC_vs_Control = "BMI",
  eoCRC_vs_loCRC   = c("BMI", "age")
)

# ── 5. Run all 6 comparisons (3 comparisons × 2 feature sets) ─────────────────
results_list <- list()

for (cmp in comparisons) {
  tag      <- cmp$tag
  covs     <- covariate_map[[tag]]
  meta_sub <- meta %>% filter(maaslin_group %in% c(cmp$case, cmp$ref))

  cat(sprintf("\n=== %s (n=%d) ===\n", tag, nrow(meta_sub)))

  # Species
  out_dir_sp <- file.path("maaslin2_results", paste0(tag, "_species"))
  fit_sp <- run_maaslin(species_prop, meta_sub, cmp$case, cmp$ref,
                        out_dir_sp, "species", covs)
  results_list[[paste0(tag, "_species")]] <- fit_sp$results

  # Pathways
  out_dir_pa <- file.path("maaslin2_results", paste0(tag, "_pathways"))
  fit_pa <- run_maaslin(path_filtered, meta_sub, cmp$case, cmp$ref,
                        out_dir_pa, "pathways", covs)
  results_list[[paste0(tag, "_pathways")]] <- fit_pa$results
}

# ── 6. Compile and save summary tables ────────────────────────────────────────
cat("\n\n=== Saving consolidated results tables ===\n")
dir.create("maaslin2_results/summaries", showWarnings = FALSE, recursive = TRUE)

all_results <- bind_rows(
  lapply(names(results_list), function(nm) {
    df <- results_list[[nm]]
    parts <- strsplit(nm, "_")[[1]]
    # tag is everything except last word (species/pathways)
    feat_type <- tail(parts, 1)
    cmp_tag   <- paste(head(parts, -1), collapse="_")
    df$comparison  <- cmp_tag
    df$feature_set <- feat_type
    df
  })
) %>%
  filter(metadata == "is_case") %>%          # keep only the group contrast row
  select(comparison, feature_set, feature, coef, stderr, pval, qval, N, N.not.zero)

write.csv(all_results,
          "maaslin2_results/summaries/all_associations.csv",
          row.names = FALSE)

# Significant hits only (q < 0.25)
sig_results <- all_results %>% filter(qval < 0.25)
write.csv(sig_results,
          "maaslin2_results/summaries/significant_q0.25.csv",
          row.names = FALSE)

cat(sprintf("  Total associations tested:  %d\n", nrow(all_results)))
cat(sprintf("  Significant (q<0.25):       %d\n", nrow(sig_results)))

# ── 7. Report: top 10 per comparison ─────────────────────────────────────────
cat("\n\n=== TOP 10 SIGNIFICANT FEATURES PER COMPARISON ===\n")
cat("(ranked by q-value, then |coef|; coef > 0 = enriched in case group)\n\n")

for (cmp in comparisons) {
  tag <- cmp$tag
  for (ftype in c("species", "pathways")) {
    key <- paste0(tag, "_", ftype)
    df  <- results_list[[key]] %>%
      filter(metadata == "is_case") %>%
      arrange(qval, desc(abs(coef))) %>%
      head(10)

    cat(sprintf("--- %s | %s (q<0.25: %d features) ---\n",
                tag, ftype,
                sum(results_list[[key]]$metadata == "is_case" &
                      results_list[[key]]$qval < 0.25, na.rm=TRUE)))

    if (nrow(df) == 0) {
      cat("  No associations at q<0.25\n\n")
    } else {
      df_print <- df %>%
        select(feature, coef, pval, qval, N.not.zero) %>%
        mutate(direction = ifelse(coef > 0,
                                  paste0("↑ in ", cmp$case),
                                  paste0("↓ in ", cmp$case)),
               coef      = round(coef, 4),
               pval      = signif(pval, 3),
               qval      = signif(qval, 3))
      print(df_print, row.names = FALSE)
      cat("\n")
    }
  }
}

# ── 8. Sample counts per group per comparison ────────────────────────────────
cat("\n=== SAMPLE COUNTS ===\n")
for (cmp in comparisons) {
  meta_sub <- meta %>% filter(maaslin_group %in% c(cmp$case, cmp$ref))
  cat(sprintf("\n%s:\n", cmp$tag))
  tbl <- table(meta_sub$maaslin_group, meta_sub$study_name)
  print(addmargins(tbl))
}

cat("\nStep 2 complete. Results saved in maaslin2_results/\n")
