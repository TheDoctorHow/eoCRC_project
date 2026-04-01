#!/usr/bin/env Rscript
# Step 2b: MaAsLin2 — corrected age-stratified grouping
#
# Groups:
#   eoCRC        : CRC, age <  50  (n=78)
#   loCRC        : CRC, age >= 50  (n=562)
#   young_ctrl   : healthy, age <  50  (n=121)
#   older_ctrl   : healthy, age >= 50  (n=521)
#
# Comparisons:
#   1. eoCRC vs young_ctrl   — fixed: group                   random: study_name
#   2. loCRC vs older_ctrl   — fixed: group + age (continuous) random: study_name
#   3. eoCRC vs loCRC        — fixed: group + age (continuous) random: study_name
#
# Feature sets: species_filtered (÷100 then CLR) and path_filtered (CLR)

suppressPackageStartupMessages({
  library(Maaslin2)
  library(dplyr)
})

setwd("/home/yugiu/eoCRC_analysis")

# ── 1. Load data ──────────────────────────────────────────────────────────────
cat("Loading data...\n")
meta             <- readRDS("all_samples_metadata.rds")
species_filtered <- readRDS("species_filtered.rds")
path_filtered    <- readRDS("path_filtered.rds")

species_prop <- species_filtered / 100   # % → proportion before CLR

stopifnot(all(rownames(species_prop)  == meta$sample_id))
stopifnot(all(rownames(path_filtered) == meta$sample_id))

# ── 2. Assign groups ──────────────────────────────────────────────────────────
meta <- meta %>%
  mutate(
    maaslin_group = case_when(
      group == "CRC"     & age <  50 ~ "eoCRC",
      group == "CRC"     & age >= 50 ~ "loCRC",
      group == "control" & age <  50 ~ "young_ctrl",
      group == "control" & age >= 50 ~ "older_ctrl"
    )
  )

cat("\nGroup counts:\n")
print(table(meta$maaslin_group))

# ── 3. Helper ─────────────────────────────────────────────────────────────────
run_maaslin <- function(feat_mat, meta_sub,
                        case_label, ref_label,
                        fixed_effects,          # character vector, NOT including study_name
                        output_dir, feat_type) {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  df_meta <- meta_sub %>%
    select(sample_id, study_name, maaslin_group, age, BMI) %>%
    mutate(group_bin = as.numeric(maaslin_group == case_label)) %>%
    as.data.frame()
  rownames(df_meta) <- df_meta$sample_id
  df_meta$sample_id <- NULL
  df_meta$maaslin_group <- NULL

  # Replace "group" placeholder in fixed_effects with actual column name
  fe <- sub("^group$", "group_bin", fixed_effects)

  feat_sub <- feat_mat[rownames(df_meta), , drop = FALSE]

  cat(sprintf(
    "  [%s | %s] %s vs %s — %d samples, %d features\n  fixed: %s | random: study_name\n",
    feat_type, basename(output_dir),
    case_label, ref_label,
    nrow(feat_sub), ncol(feat_sub),
    paste(fe, collapse = " + ")
  ))

  fit <- Maaslin2(
    input_data       = feat_sub,
    input_metadata   = df_meta,
    output           = output_dir,
    fixed_effects    = fe,
    random_effects   = "study_name",
    normalization    = "CLR",
    transform        = "NONE",
    analysis_method  = "LM",
    max_significance = 0.25,
    correction       = "BH",
    min_prevalence   = 0,
    min_abundance    = 0,
    cores            = 4,
    plot_heatmap     = FALSE,
    plot_scatter     = FALSE
  )
  invisible(fit)
}

# ── 4. Comparison specs ───────────────────────────────────────────────────────
comparisons <- list(
  list(
    tag    = "eoCRC_vs_YoungCtrl",
    case   = "eoCRC",
    ref    = "young_ctrl",
    fixed  = "group"                  # no age covariate (both < 50)
  ),
  list(
    tag    = "loCRC_vs_OlderCtrl",
    case   = "loCRC",
    ref    = "older_ctrl",
    fixed  = c("group", "age")        # age continuous within older stratum
  ),
  list(
    tag    = "eoCRC_vs_loCRC",
    case   = "eoCRC",
    ref    = "loCRC",
    fixed  = c("group", "age")        # age essential to separate the groups
  )
)

# ── 5. Run all 6 MaAsLin2 models ─────────────────────────────────────────────
results_list <- list()
base_dir <- "maaslin2_results_v2"

for (cmp in comparisons) {
  meta_sub <- meta %>% filter(maaslin_group %in% c(cmp$case, cmp$ref))
  cat(sprintf("\n=== %s (n=%d) ===\n", cmp$tag, nrow(meta_sub)))

  for (ftype in c("species", "pathways")) {
    feat_mat <- if (ftype == "species") species_prop else path_filtered
    out_dir  <- file.path(base_dir, paste0(cmp$tag, "_", ftype))
    key      <- paste0(cmp$tag, "_", ftype)

    fit <- run_maaslin(feat_mat, meta_sub,
                       cmp$case, cmp$ref,
                       cmp$fixed, out_dir, ftype)
    results_list[[key]] <- fit$results
  }
}

# ── 6. Consolidate and save ───────────────────────────────────────────────────
cat("\n\n=== Consolidating results ===\n")
dir.create(file.path(base_dir, "summaries"), showWarnings = FALSE, recursive = TRUE)

all_results <- bind_rows(
  lapply(names(results_list), function(nm) {
    df    <- results_list[[nm]]
    parts <- strsplit(nm, "_")[[1]]
    ftype <- tail(parts, 1)
    ctag  <- paste(head(parts, -1), collapse = "_")
    df$comparison  <- ctag
    df$feature_set <- ftype
    df
  })
) %>%
  filter(metadata == "group_bin") %>%
  select(comparison, feature_set, feature, coef, stderr, pval, qval, N, N.not.zero)

write.csv(all_results,
          file.path(base_dir, "summaries", "all_associations.csv"),
          row.names = FALSE)

sig_results <- all_results %>% filter(qval < 0.25)
write.csv(sig_results,
          file.path(base_dir, "summaries", "significant_q0.25.csv"),
          row.names = FALSE)

cat(sprintf("  Total associations: %d\n", nrow(all_results)))
cat(sprintf("  Significant q<0.25: %d\n", nrow(sig_results)))
cat(sprintf("  Saved: %s/summaries/\n", base_dir))

# ── 7. Top-10 report ──────────────────────────────────────────────────────────
cat("\n\n=== TOP 10 FEATURES PER COMPARISON (ranked by q-value) ===\n")
cat("  coef > 0 = enriched in case group (first named)\n\n")

for (cmp in comparisons) {
  for (ftype in c("species", "pathways")) {
    key   <- paste0(cmp$tag, "_", ftype)
    df    <- results_list[[key]] %>%
      filter(metadata == "group_bin") %>%
      arrange(qval, desc(abs(coef)))
    n_sig <- sum(df$qval < 0.25, na.rm = TRUE)

    cat(sprintf("--- %s | %s  (q<0.25: %d / %d tested) ---\n",
                cmp$tag, ftype, n_sig, nrow(df)))

    top10 <- head(df, 10) %>%
      mutate(
        direction = ifelse(coef > 0,
                           paste0("↑ ", cmp$case),
                           paste0("↓ ", cmp$case)),
        coef      = round(coef, 5),
        pval      = signif(pval, 3),
        qval      = signif(qval, 3)
      ) %>%
      select(feature, coef, pval, qval, N.not.zero, direction)

    print(top10, row.names = FALSE)
    cat("\n")
  }
}

# ── 8. Sample counts per comparison ──────────────────────────────────────────
cat("\n=== SAMPLE COUNTS PER COMPARISON ===\n")
for (cmp in comparisons) {
  meta_sub <- meta %>% filter(maaslin_group %in% c(cmp$case, cmp$ref))
  cat(sprintf("\n%s:\n", cmp$tag))
  tbl <- table(meta_sub$maaslin_group, meta_sub$study_name)
  print(addmargins(tbl))
}

cat("\nStep 2b complete.\n")
