#!/usr/bin/env Rscript
# Step 1: Extract and prevalence-filter species and pathway matrices
# Output: filtered matrices as RDS + CSV, plus a summary report

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(TreeSummarizedExperiment)
})

setwd("/home/yugiu/eoCRC_analysis")

cat("=== Step 1: Data Extraction & Prevalence Filtering ===\n\n")

# ── Load objects ──────────────────────────────────────────────────────────────
cat("Loading RDS objects...\n")
tse_full  <- readRDS("tse_relabund_full.rds")
tse_path  <- readRDS("tse_pathways_full.rds")
meta      <- readRDS("all_samples_metadata.rds")

cat("  tse_relabund_full dim:", paste(dim(tse_full),  collapse=" x "), "\n")
cat("  tse_pathways_full dim:", paste(dim(tse_path),  collapse=" x "), "\n")
cat("  metadata rows:        ", nrow(meta), "\n\n")

# ── Inspect metadata ──────────────────────────────────────────────────────────
cat("Metadata columns:", paste(colnames(meta), collapse=", "), "\n")
cat("Group distribution:\n")
print(table(meta$group, meta$age_group, useNA="ifany"))
cat("\nStudy counts:\n")
print(table(meta$study_name))
cat("\n")

# ── Extract species matrix (samples x species) ────────────────────────────────
cat("Extracting species matrix...\n")
species_mat_raw <- t(assay(tse_full, "relative_abundance"))
cat("  Raw species matrix:", paste(dim(species_mat_raw), collapse=" x "),
    "(samples x species)\n")

# ── Extract pathway matrix — unstratified only ────────────────────────────────
cat("Extracting pathway matrix...\n")
path_mat_raw <- t(assay(tse_path, "pathway_abundance"))
cat("  Raw pathway matrix (before filter):", paste(dim(path_mat_raw), collapse=" x "), "\n")

# Remove stratified (pipes), UNMAPPED, UNINTEGRATED
keep_paths <- !grepl("\\|", colnames(path_mat_raw)) &
              !grepl("UNMAPPED|UNINTEGRATED", colnames(path_mat_raw))
path_mat <- path_mat_raw[, keep_paths]
cat("  After removing stratified/UNMAPPED/UNINTEGRATED:",
    paste(dim(path_mat), collapse=" x "), "\n\n")

# ── Prevalence filtering (≥10% of samples) ───────────────────────────────────
prev_filter <- function(mat, min_prev = 0.10, label = "") {
  prevalence  <- colMeans(mat > 0)
  keep        <- prevalence >= min_prev
  cat(sprintf("  %s: %d features before → %d after prevalence filter (%.0f%% kept)\n",
              label, ncol(mat), sum(keep), 100 * mean(keep)))
  list(filtered = mat[, keep], prevalence = prevalence, kept = keep)
}

cat("Applying prevalence filter (min 10% of samples)...\n")
sp_res   <- prev_filter(species_mat_raw, label = "Species")
path_res <- prev_filter(path_mat,        label = "Pathways")

species_filtered <- sp_res$filtered
path_filtered    <- path_res$filtered

cat("\nFinal dimensions:\n")
cat("  species_filtered:", paste(dim(species_filtered), collapse=" x "),
    "(samples x species)\n")
cat("  path_filtered:   ", paste(dim(path_filtered),    collapse=" x "),
    "(samples x pathways)\n\n")

# ── Quick sanity checks ───────────────────────────────────────────────────────
# Confirm row order matches metadata
stopifnot(all(rownames(species_filtered) == meta$sample_id |
              nrow(species_filtered) == nrow(meta)))
cat("Row order check passed (", nrow(species_filtered), "samples )\n\n")

# Relative abundance sums (should be ≈1 per sample)
row_sums <- rowSums(species_mat_raw)
cat("Relative abundance row sums — min:", round(min(row_sums), 4),
    " median:", round(median(row_sums), 4),
    " max:", round(max(row_sums), 4), "\n\n")

# ── Save filtered matrices ────────────────────────────────────────────────────
cat("Saving filtered matrices...\n")
saveRDS(species_filtered, "species_filtered.rds")
saveRDS(path_filtered,    "path_filtered.rds")

# Save as CSV too (large, but useful for inspection)
write.csv(species_filtered, "species_filtered.csv")
write.csv(path_filtered,    "path_filtered.csv")
cat("  Saved: species_filtered.rds / .csv\n")
cat("  Saved: path_filtered.rds / .csv\n\n")

# ── Save prevalence vectors ───────────────────────────────────────────────────
saveRDS(sp_res$prevalence,   "species_prevalence.rds")
saveRDS(path_res$prevalence, "path_prevalence.rds")
cat("  Saved: species_prevalence.rds, path_prevalence.rds\n\n")

# ── Summary report ────────────────────────────────────────────────────────────
cat("=== SUMMARY ===\n")
cat(sprintf("Species matrix:  %d samples x %d species (raw) → %d species (filtered)\n",
            nrow(species_mat_raw), ncol(species_mat_raw), ncol(species_filtered)))
cat(sprintf("Pathway matrix:  %d samples x %d pathways (unstratified) → %d pathways (filtered)\n",
            nrow(path_mat), ncol(path_mat), ncol(path_filtered)))

# Key taxa check
key_taxa <- c("Fusobacterium_nucleatum", "Peptostreptococcus_anaerobius",
              "Escherichia_coli", "Faecalibacterium_prausnitzii")
cat("\nKey taxa in filtered matrix:\n")
for (tx in key_taxa) {
  match <- grep(tx, colnames(species_filtered), value=TRUE)
  if (length(match) > 0) {
    cat(sprintf("  ✓ %s → %s\n", tx, paste(match, collapse=", ")))
  } else {
    cat(sprintf("  ✗ %s — NOT found (may be filtered out or name differs)\n", tx))
  }
}

cat("\nStep 1 complete.\n")
