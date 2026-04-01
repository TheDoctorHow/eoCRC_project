#!/usr/bin/env Rscript
# Step 2c: MaAsLin2 v3 — proper CLR pre-computation
#
# Fix: MaAsLin2 v1.24.1 uses chemometrics::clr(x + 1) internally, which for
#      proportions (x << 1) reduces to arithmetic mean-centering, NOT true CLR.
#      Solution: pre-compute CLR externally with compositions::clr(x + 1e-6),
#      then pass to MaAsLin2 with normalization='NONE'.
#
# Comparisons (same as v2):
#   1. eoCRC vs young_ctrl
#   2. loCRC vs older_ctrl
#   3. eoCRC vs loCRC
#
# Feature sets: species (proportions / 100) and pathways (CPM)

suppressPackageStartupMessages({
  library(Maaslin2)
  library(compositions)
  library(dplyr)
})

setwd("/home/yugiu/eoCRC_analysis")

PSEUDO <- 1e-6
BASE   <- "maaslin2_results_v3"
dir.create(file.path(BASE, "summaries"), showWarnings = FALSE, recursive = TRUE)

# ── 1. Load data ───────────────────────────────────────────────────────────────
cat("Loading data...\n")
meta             <- readRDS("all_samples_metadata.rds")
species_filtered <- readRDS("species_filtered.rds")   # proportions in %, 0-100
path_filtered    <- readRDS("path_filtered.rds")       # pathway CPM

# Convert species to proportions (0-1) before CLR
species_prop <- species_filtered / 100

meta <- meta %>%
  mutate(maaslin_group = case_when(
    group == "CRC"     & age <  50 ~ "eoCRC",
    group == "CRC"     & age >= 50 ~ "loCRC",
    group == "control" & age <  50 ~ "young_ctrl",
    group == "control" & age >= 50 ~ "older_ctrl"
  ))

cat("Group counts:\n"); print(table(meta$maaslin_group))

# ── 2. CLR helper ──────────────────────────────────────────────────────────────
# Apply CLR to a sample SUBSET so geometric mean is computed within-comparison.
# compositions::clr(x) computes log(x_j) - mean(log(x_j)) row-wise.
clr_subset <- function(mat, pseudo = PSEUDO) {
  as.matrix(compositions::clr(mat + pseudo))
}

# ── 3. MaAsLin2 runner ────────────────────────────────────────────────────────
run_maaslin_v3 <- function(clr_mat, meta_sub, case_label,
                           fixed_effects, output_dir, feat_type) {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  df_meta <- meta_sub %>%
    select(sample_id, study_name, maaslin_group, age) %>%
    mutate(group_bin = as.numeric(maaslin_group == case_label)) %>%
    as.data.frame()
  rownames(df_meta) <- df_meta$sample_id
  df_meta$sample_id     <- NULL
  df_meta$maaslin_group <- NULL

  fe <- sub("^group$", "group_bin", fixed_effects)

  feat_sub <- clr_mat[rownames(df_meta), , drop = FALSE]

  cat(sprintf(
    "  [%s] %d samples, %d features | fixed: %s | random: study_name\n",
    basename(output_dir), nrow(feat_sub), ncol(feat_sub),
    paste(fe, collapse = " + ")
  ))

  Maaslin2(
    input_data       = feat_sub,
    input_metadata   = df_meta,
    output           = output_dir,
    fixed_effects    = fe,
    random_effects   = "study_name",
    normalization    = "NONE",      # CLR already applied externally
    transform        = "NONE",
    analysis_method  = "LM",
    standardize      = TRUE,        # standardise continuous metadata (group_bin, age)
    max_significance = 0.25,
    correction       = "BH",
    min_prevalence   = 0,
    min_abundance    = 0,
    cores            = 4,
    plot_heatmap     = FALSE,
    plot_scatter     = FALSE
  )
}

# ── 4. Comparison specs ────────────────────────────────────────────────────────
comparisons <- list(
  list(tag="eoCRC_vs_YoungCtrl", case="eoCRC",   ref="young_ctrl", fixed="group"),
  list(tag="loCRC_vs_OlderCtrl", case="loCRC",   ref="older_ctrl", fixed=c("group","age")),
  list(tag="eoCRC_vs_loCRC",     case="eoCRC",   ref="loCRC",      fixed=c("group","age"))
)

# ── 5. Run all 6 models ────────────────────────────────────────────────────────
results_list <- list()

for (cmp in comparisons) {
  meta_sub <- meta %>% filter(maaslin_group %in% c(cmp$case, cmp$ref))
  cat(sprintf("\n=== %s (n=%d) ===\n", cmp$tag, nrow(meta_sub)))

  for (ftype in c("species", "pathways")) {
    raw_mat <- if (ftype == "species") species_prop else path_filtered
    # CLR computed on this comparison's samples only
    clr_mat <- clr_subset(raw_mat[meta_sub$sample_id, , drop = FALSE])

    cat(sprintf("  CLR range [%s]: %.3f to %.3f\n",
                ftype, min(clr_mat), max(clr_mat)))

    out_dir <- file.path(BASE, paste0(cmp$tag, "_", ftype))
    key     <- paste0(cmp$tag, "_", ftype)

    fit <- run_maaslin_v3(clr_mat, meta_sub, cmp$case, cmp$fixed, out_dir, ftype)
    results_list[[key]] <- fit$results
  }
}

# ── 6. Consolidate — read from written TSV files (column names are stable) ─────
cat("\n\n=== Consolidating results ===\n")

all_results <- bind_rows(lapply(names(results_list), function(nm) {
  parts   <- strsplit(nm, "_")[[1]]
  ftype   <- tail(parts, 1)
  ctag    <- paste(head(parts, -1), collapse = "_")
  tsv     <- file.path(BASE, nm, "all_results.tsv")
  df      <- read.delim(tsv, stringsAsFactors = FALSE)
  df$comparison  <- ctag
  df$feature_set <- ftype
  df
})) %>%
  filter(metadata == "group_bin") %>%
  select(comparison, feature_set, feature, coef, stderr, pval, qval, N, N.not.0)

write.csv(all_results,
          file.path(BASE, "summaries", "all_associations.csv"),
          row.names = FALSE)

sig_results <- all_results %>% filter(qval < 0.25)
write.csv(sig_results,
          file.path(BASE, "summaries", "significant_q0.25.csv"),
          row.names = FALSE)

cat(sprintf("  Total associations: %d\n", nrow(all_results)))
cat(sprintf("  Significant q<0.25: %d\n", nrow(sig_results)))

# ── 7. Consistency check vs v2 ────────────────────────────────────────────────
cat("\n\n=== CONSISTENCY CHECK: v3 vs v2 (eoCRC_vs_YoungCtrl, species) ===\n")

v2 <- read.csv("maaslin2_results_v2/summaries/significant_q0.25.csv",
               stringsAsFactors = FALSE) %>%
  filter(comparison == "eoCRC_vs_YoungCtrl", feature_set == "species") %>%
  select(feature, coef_v2 = coef, qval_v2 = qval)

v3 <- sig_results %>%
  filter(comparison == "eoCRC_vs_YoungCtrl", feature_set == "species") %>%
  select(feature, coef_v3 = coef, qval_v3 = qval)

check <- full_join(v2, v3, by = "feature") %>%
  mutate(
    same_direction = sign(coef_v2) == sign(coef_v3),
    both_sig       = !is.na(coef_v2) & !is.na(coef_v3)
  )

cat("Species significant in v2:\n")
cat(paste(" ", v2$feature, collapse = "\n"), "\n\n")
cat("Species significant in v3:\n")
cat(paste(" ", v3$feature, collapse = "\n"), "\n\n")
cat("Direction consistent for shared hits:", all(check$same_direction, na.rm=TRUE), "\n")

# ── 8. Top 10 report for primary comparison ───────────────────────────────────
cat("\n\n=== TOP 10: eoCRC_vs_YoungCtrl | species (v3, proper CLR) ===\n")
cat("coef = CLR effect size per SD(group_bin)\n")
cat("  [multiply by ~2.05 to get absolute CLR difference eoCRC - young_ctrl]\n\n")

top10 <- read.delim(file.path(BASE, "eoCRC_vs_YoungCtrl_species", "all_results.tsv"),
                    stringsAsFactors = FALSE) %>%
  filter(metadata == "group_bin") %>%
  arrange(qval, desc(abs(coef))) %>%
  head(10) %>%
  mutate(
    species   = gsub("^species\\.", "", feature),
    species   = gsub("\\.", " ", species),
    direction = ifelse(coef > 0, paste0("UP in eoCRC"), paste0("DOWN in eoCRC")),
    coef      = round(coef, 4),
    stderr    = round(stderr, 4),
    pval      = signif(pval, 3),
    qval      = signif(qval, 3)
  ) %>%
  select(species, coef, stderr, pval, qval, direction)

print(top10, row.names = FALSE)

cat("\n\nStep 2c complete. Results in maaslin2_results_v3/\n")
