#!/usr/bin/env Rscript
# Step 4: SHAP interpretation of final Combined RF model
# Model: RF (ntree=500) trained on ALL eoCRC + young controls (n=199)
# Features: combined species + pathways, CLR within-sample, top 50 by variance
# Backend: fastshap (approximate Shapley values, nsim=100)

suppressPackageStartupMessages({
  library(randomForest)
  library(fastshap)
  library(shapviz)
  library(compositions)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

setwd("/home/yugiu/eoCRC_analysis")
set.seed(42)
dir.create("figures", showWarnings = FALSE)

N_TOP    <- 50
N_SHAP   <- 20   # top features to show in plots
NSIM     <- 100  # fastshap simulations (accuracy vs speed)
NTREE    <- 500

# ── 1. Load and prepare data ──────────────────────────────────────────────────
cat("Preparing data...\n")
meta <- readRDS("all_samples_metadata.rds") %>%
  mutate(maaslin_group = case_when(
    group == "CRC"     & age <  50 ~ "eoCRC",
    group == "CRC"     & age >= 50 ~ "loCRC",
    group == "control" & age <  50 ~ "young_ctrl",
    group == "control" & age >= 50 ~ "older_ctrl"
  ))

sp_mat <- readRDS("species_filtered.rds")
pa_mat <- readRDS("path_filtered.rds")

idx    <- meta$maaslin_group %in% c("eoCRC", "young_ctrl")
y      <- as.integer(meta$maaslin_group[idx] == "eoCRC")  # 1=eoCRC
sp_sub <- sp_mat[idx, ];  pa_sub <- pa_mat[idx, ]

cat(sprintf("  n=%d  eoCRC=%d  young_ctrl=%d\n\n", sum(idx), sum(y), sum(y==0)))

# ── 2. CLR + combined feature matrix (same as LOSO-CV but on full data) ───────
cat("Applying CLR transformation and selecting top-50 features...\n")
pseudo <- 1e-6
sp_clr <- as.matrix(compositions::clr(sp_sub + pseudo))
pa_clr <- as.matrix(compositions::clr(pa_sub + pseudo))
X_comb <- cbind(sp_clr, pa_clr)     # 199 x 639

# Variance-based feature selection on full data
vars    <- apply(X_comb, 2, var)
top50   <- names(sort(vars, decreasing = TRUE))[seq_len(N_TOP)]
X_final <- X_comb[, top50]

cat(sprintf("  Combined matrix: %d x %d -> selected top %d by variance\n",
            nrow(X_comb), ncol(X_comb), N_TOP))
cat("  Top 5 features:", paste(head(top50, 5), collapse=", "), "\n\n")

# ── 3. Train final RF on all data ─────────────────────────────────────────────
cat(sprintf("Training final RF (ntree=%d) on all %d samples...\n", NTREE, nrow(X_final)))
final_rf <- randomForest(x = X_final, y = as.factor(y), ntree = NTREE, importance = TRUE)
cat(sprintf("  OOB error: %.1f%%\n\n", final_rf$err.rate[NTREE, "OOB"] * 100))

# ── 4. Compute SHAP values via fastshap ───────────────────────────────────────
cat(sprintf("Computing SHAP values (fastshap, nsim=%d)...\n", NSIM))
cat("  This takes ~2-3 minutes...\n")

pred_fn <- function(model, newdata)
  predict(model, newdata, type = "prob")[, "1"]

shap_raw <- fastshap::explain(
  final_rf,
  X            = as.data.frame(X_final),
  pred_wrapper = pred_fn,
  nsim         = NSIM,
  adjust       = TRUE
)
# shapviz object
shap_obj <- shapviz(shap_raw, X = X_final)

# Save SHAP matrix
shap_mat <- as.data.frame(shap_raw)
rownames(shap_mat) <- rownames(X_final)
saveRDS(shap_mat, "shap_values.rds")
write.csv(shap_mat, "shap_values.csv")
cat(sprintf("  SHAP matrix saved: %d samples x %d features\n\n", nrow(shap_mat), ncol(shap_mat)))

# ── 5. Identify top features ──────────────────────────────────────────────────
mean_abs_shap <- colMeans(abs(shap_mat))
top20_names   <- names(sort(mean_abs_shap, decreasing = TRUE))[seq_len(N_SHAP)]

# Feature type: species or pathway?
feat_type <- ifelse(grepl("^species:", top20_names), "species", "pathway")

# Direction: mean SHAP > 0 → enriched in eoCRC
mean_shap  <- colMeans(shap_mat)
direction  <- ifelse(mean_shap[top20_names] > 0, "↑ eoCRC", "↓ eoCRC")

# ── 6. Cross-reference with MaAsLin2 ─────────────────────────────────────────
cat("Cross-referencing with MaAsLin2 results...\n")
maaslin_sig <- read.csv("maaslin2_results_v2/summaries/significant_q0.25.csv",
                        stringsAsFactors = FALSE) %>%
  filter(comparison == "eoCRC_vs_YoungCtrl")

# Normalize names for comparison: replace non-alphanumeric with "_", lowercase
norm <- function(x) tolower(gsub("[^A-Za-z0-9]", "_", x))

shap_norm    <- norm(top20_names)
maaslin_norm <- norm(maaslin_sig$feature)

in_maaslin <- shap_norm %in% maaslin_norm

# Build annotation for matched MaAsLin2 hits
matched_q <- sapply(seq_along(top20_names), function(i) {
  if (!in_maaslin[i]) return(NA_real_)
  idx_m <- which(maaslin_norm == shap_norm[i])[1]
  maaslin_sig$qval[idx_m]
})
matched_coef <- sapply(seq_along(top20_names), function(i) {
  if (!in_maaslin[i]) return(NA_real_)
  idx_m <- which(maaslin_norm == shap_norm[i])[1]
  maaslin_sig$coef[idx_m]
})

# ── 7. Summary table ──────────────────────────────────────────────────────────
top20_df <- data.frame(
  rank         = seq_len(N_SHAP),
  feature      = top20_names,
  feat_type    = feat_type,
  mean_abs_SHAP= round(mean_abs_shap[top20_names], 5),
  direction    = direction,
  in_maaslin2  = ifelse(in_maaslin, "YES", "no"),
  maaslin2_q   = round(matched_q,   4),
  maaslin2_coef= round(matched_coef, 5),
  validated    = ifelse(in_maaslin, "DOUBLY VALIDATED", "SHAP only"),
  stringsAsFactors = FALSE
)

# Short display name (strip "species:" prefix, truncate pathways)
top20_df$short_name <- gsub("^species:", "", top20_df$feature)
top20_df$short_name <- gsub("^.*PWY[._-]?([0-9]+).*$", "PWY-\\1", top20_df$short_name)
top20_df$short_name <- substr(top20_df$short_name, 1, 40)

saveRDS(top20_df, "shap_top20_summary.rds")
write.csv(top20_df, "shap_top20_summary.csv", row.names = FALSE)

cat("\n=== TOP 20 SHAP FEATURES ===\n")
print(top20_df[, c("rank","feature","feat_type","mean_abs_SHAP",
                   "direction","in_maaslin2","maaslin2_q","validated")],
      row.names = FALSE)

cat(sprintf("\nDoubly validated (in both SHAP top-20 AND MaAsLin2 q<0.25): %d\n",
            sum(in_maaslin)))
cat(sprintf("Species in top 20: %d  |  Pathways: %d\n\n",
            sum(feat_type == "species"), sum(feat_type == "pathway")))

# ── 8. Figures ────────────────────────────────────────────────────────────────
cat("Generating figures...\n")

# Helper: clean feature labels for axis
clean_label <- function(x) {
  x <- gsub("^species:", "", x)
  x <- gsub("\\.", " ", x)
  x <- gsub("_", " ", x)
  # Abbreviate long pathway names
  ifelse(nchar(x) > 45, paste0(substr(x, 1, 43), "…"), x)
}

# Restrict shapviz object to top-20 features for cleaner plots
shap_top20 <- shap_obj[, top20_names]

# Rename features to clean labels in the shapviz object
colnames(shap_top20$S) <- clean_label(top20_names)
colnames(shap_top20$X) <- clean_label(top20_names)

# ── Fig 1: Beeswarm plot (top 20, high=red, low=blue) ────────────────────────
fig_beeswarm <- sv_importance(
  shap_top20,
  kind       = "beeswarm",
  max_display= N_SHAP,
  color_bar_title = "Feature value\n(CLR)"
) +
  labs(
    title    = "SHAP beeswarm — eoCRC vs Young Controls",
    subtitle = sprintf("Combined RF (ntree=%d) trained on all n=199 samples | top 20 features", NTREE),
    x        = "SHAP value  (positive = ↑ eoCRC probability)"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12))

ggsave("figures/fig_shap_beeswarm.png", fig_beeswarm,
       width = 9, height = 7, dpi = 300)
cat("  Saved: figures/fig_shap_beeswarm.png\n")

# ── Fig 2: Bar plot (mean |SHAP|) ─────────────────────────────────────────────
# Build manually for full control + doubly-validated annotation
bar_df <- top20_df %>%
  mutate(
    label     = clean_label(feature),
    label     = factor(label, levels = rev(clean_label(top20_names))),
    fill_col  = case_when(
      validated == "DOUBLY VALIDATED" & feat_type == "species"  ~ "Species (validated)",
      validated == "DOUBLY VALIDATED" & feat_type == "pathway"  ~ "Pathway (validated)",
      feat_type == "species"  ~ "Species",
      TRUE                    ~ "Pathway"
    )
  )

fill_palette <- c(
  "Species (validated)" = "#C0392B",
  "Pathway (validated)" = "#1A5276",
  "Species"             = "#F1948A",
  "Pathway"             = "#7FB3D3"
)

fig_bar <- ggplot(bar_df, aes(x = mean_abs_SHAP, y = label, fill = fill_col)) +
  geom_col() +
  geom_text(aes(label = ifelse(validated == "DOUBLY VALIDATED", "★", "")),
            hjust = -0.2, size = 4, color = "black") +
  scale_fill_manual(values = fill_palette, name = NULL) +
  labs(
    title    = "Mean |SHAP| importance — top 20 features",
    subtitle = "★ = doubly validated (SHAP top-20 + MaAsLin2 q<0.25)",
    x        = "Mean |SHAP value|", y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title   = element_text(face = "bold", size = 12),
    legend.position = "bottom"
  )

ggsave("figures/fig_shap_bar.png", fig_bar,
       width = 9, height = 7, dpi = 300)
cat("  Saved: figures/fig_shap_bar.png\n")

# ── Fig 3: Dependence plots — top 5 features ─────────────────────────────────
top5 <- top20_names[1:5]

dep_plots <- lapply(seq_along(top5), function(i) {
  fn <- top5[i]
  cl <- clean_label(fn)
  # Use full shap_obj (not top20 subset) for interaction colouring
  sv_dependence(shap_obj, v = fn, color_var = "auto") +
    labs(
      title    = cl,
      x        = paste0(cl, "\n(CLR-transformed abundance)"),
      y        = "SHAP value"
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10),
          legend.key.size = unit(0.4, "cm"))
})

fig_dep <- wrap_plots(dep_plots, ncol = 2) +
  plot_annotation(
    title    = "SHAP dependence plots — top 5 features",
    subtitle = "SHAP value vs CLR-transformed abundance | color = most interacting feature",
    theme    = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave("figures/fig_shap_dependence.png", fig_dep,
       width = 12, height = 10, dpi = 300)
cat("  Saved: figures/fig_shap_dependence.png\n")

# ── Fig 4: Combined beeswarm + bar (presentation-ready) ───────────────────────
fig_combined <- fig_beeswarm + fig_bar +
  plot_layout(widths = c(1.2, 1)) +
  plot_annotation(
    title  = "SHAP feature importance — eoCRC microbiome signature",
    theme  = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave("figures/fig_shap_combined.png", fig_combined,
       width = 16, height = 7, dpi = 300)
cat("  Saved: figures/fig_shap_combined.png\n\n")

# ── 9. Literature cross-reference for top 5 ───────────────────────────────────
lit_notes <- list(
  "species:Parvimonas micra"    = "Oral-origin pathobiont; consistently enriched in CRC across studies; associated with left-sided and rectal tumors; possible eoCRC-specific marker (Komiya et al. 2019)",
  "species:Blautia obeum"       = "Ambiguous role; some studies show enrichment in CRC; may reflect dysbiosis rather than direct oncogenic activity",
  "species:Gemella morbillorum" = "Oral-origin commensal; enriched in CRC in multiple cohorts; co-occurs with Fusobacterium",
  "species:Fusobacterium nucleatum" = "Most replicated CRC microbiome marker; FadA adhesin promotes Wnt signaling; associated with MSI-H CRC; more prevalent in loCRC in our data",
  "species:Faecalibacterium prausnitzii" = "Major butyrate producer; depleted in IBD and CRC; anti-inflammatory via HDAC inhibition; depletion signals loss of mucosal protection",
  "species:Bacteroides caccae"  = "Enriched in multiple CRC cohorts; metabolizes complex polysaccharides; produces pro-inflammatory LPS",
  "species:Peptostreptococcus stomatis" = "Oral-origin; strong CRC marker in Thomas et al. 2019; possible biofilm participant",
  "species:Alistipes finegoldii" = "Hydrogen sulfide producer; enriched in CRC; sulfide damages colonocyte mitochondria",
  "species:Alistipes indistinctus" = "Less characterized; same Alistipes genus as finegoldii; sulfide producer family",
  "species:[Clostridium] symbiosum" = "Associated with CRC in multiple cohorts; produces oncometabolites; butyrate-consuming competitor"
)

cat("=== LITERATURE NOTES FOR TOP 5 SHAP FEATURES ===\n\n")
for (i in 1:5) {
  fn <- top5[i]
  cat(sprintf("Rank %d: %s\n", i, fn))
  cat(sprintf("  Direction:     %s\n", direction[i]))
  cat(sprintf("  Mean |SHAP|:   %.5f\n", mean_abs_shap[fn]))
  cat(sprintf("  In MaAsLin2:   %s", ifelse(in_maaslin[i], "YES", "no")))
  if (in_maaslin[i]) cat(sprintf("  (q=%.4f, coef=%.4f)", matched_q[i], matched_coef[i]))
  cat("\n")
  note <- lit_notes[[fn]]
  if (!is.null(note)) cat(sprintf("  Literature:    %s\n", note))
  cat("\n")
}

# ── 10. Proportion species vs pathways ────────────────────────────────────────
cat("=== FEATURE TYPE BREAKDOWN ===\n")
type_tab <- table(top20_df$feat_type)
cat(sprintf("Species in top 20:  %d / %d (%.0f%%)\n",
            type_tab["species"], N_SHAP, 100*type_tab["species"]/N_SHAP))
cat(sprintf("Pathways in top 20: %d / %d (%.0f%%)\n",
            type_tab["pathway"], N_SHAP, 100*type_tab["pathway"]/N_SHAP))

validated_by_type <- table(top20_df$feat_type[in_maaslin])
cat(sprintf("\nDoubly validated species:  %d\n",
            ifelse("species" %in% names(validated_by_type), validated_by_type["species"], 0)))
cat(sprintf("Doubly validated pathways: %d\n",
            ifelse("pathway" %in% names(validated_by_type), validated_by_type["pathway"], 0)))

cat("\nStep 4 complete. Figures saved to figures/\n")
