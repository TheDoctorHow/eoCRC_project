#!/usr/bin/env Rscript
# Regenerate all 7 SHAP supplementary figures with v3 MaAsLin2 doubly-validated list
# 13 doubly-validated features (up from 8 in v2 cross-reference)
# Label coloring: enriched (red), depleted (blue), non-validated (black)
# Reference style: fig4_shap_beeswarm.png (regen_fig4.R)

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
  library(patchwork); library(ggrepel)
})

setwd("/home/yugiu/eoCRC_analysis")
dir.create("figures", showWarnings = FALSE)

# ── Palette ────────────────────────────────────────────────────────────────────
COL_UP  <- "#C0392B"   # enriched in eoCRC (red)
COL_DN  <- "#1A5276"   # depleted in eoCRC (blue)
COL_NS  <- "grey20"    # non-validated (dark grey)

theme_pres <- function(base = 12)
  theme_bw(base_size = base) +
  theme(plot.title    = element_text(face = "bold", size = base + 1),
        plot.subtitle = element_text(size = base - 1, color = "grey40"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", color = NA),
        strip.text = element_text(face = "bold"))

# ── v3 doubly-validated features (13 total) ───────────────────────────────────
V3_VALIDATED <- c(
  "species:Peptostreptococcus stomatis",
  "species:Gemella morbillorum",
  "species:Parvimonas micra",
  "species:Intestinimonas butyriciproducens",
  "species:Dialister pneumosintes",
  "species:Faecalibacterium prausnitzii",
  "species:Lactobacillus rogosae",
  "species:Alistipes finegoldii",
  "species:Anaerotruncus colihominis",
  "species:Bacteroides cellulosilyticus",
  "species:Bacteroides caccae",
  "species:Turicimonas muris",
  "species:[Clostridium] symbiosum"
)

# Depleted in eoCRC (MaAsLin2 v3 coef < 0) → blue labels
V3_DEPLETED <- c(
  "species:Faecalibacterium prausnitzii",    # coef = -0.477, q = 0.082
  "species:Lactobacillus rogosae",            # coef = -0.394, q = 0.078
  "species:Turicimonas muris"                # coef = -0.513, q = 0.024
)

# ── Helper functions ────────────────────────────────────────────────────────────
clean_sp <- function(x) {
  x <- gsub("^species:", "", x)
  x <- gsub("[._]", " ", x)
  x <- gsub("sp  ", "sp. ", x)     # restore "sp. X" notation
  x <- gsub("CAG  ", "CAG:", x)    # restore "CAG:NNN"
  ifelse(nchar(x) > 40, paste0(substr(x, 1, 38), "\u2026"), x)
}

label_for <- function(feat) {
  cl <- clean_sp(feat)
  ifelse(feat %in% V3_VALIDATED, paste0(cl, "  \u2605"), cl)
}

label_color <- function(lbl, feat_vec) {
  # lbl and feat_vec are parallel vectors (same length, same order)
  sapply(seq_along(lbl), function(i) {
    f <- feat_vec[i]
    if (!grepl("\u2605", lbl[i])) return(COL_NS)
    if (f %in% V3_DEPLETED)       return(COL_DN)
    return(COL_UP)
  })
}

# ── Load data ──────────────────────────────────────────────────────────────────
cat("Loading SHAP extended results...\n")
shap_ext   <- readRDS("shap_extended_results.rds")
shap_df    <- shap_ext$df       # top-20 summary data frame
shap_mat   <- shap_ext$shap_mat # 199 x 200 SHAP matrix
X_shap     <- shap_ext$Xsel     # 199 x 200 feature matrix

top20_feat <- shap_df$feature   # ordered by mean |SHAP|

# Update validation status in shap_df using v3 list
shap_df <- shap_df %>%
  mutate(
    validated  = ifelse(feature %in% V3_VALIDATED, "DOUBLY VALIDATED", "SHAP only"),
    depleted   = feature %in% V3_DEPLETED,
    enriched   = (feature %in% V3_VALIDATED) & !(feature %in% V3_DEPLETED),
    clean_name = clean_sp(feature),
    feat_label = label_for(feature)
  )

n_val <- sum(shap_df$validated == "DOUBLY VALIDATED")
cat(sprintf("  Doubly validated in top-20: %d / 20\n", n_val))
cat(sprintf("  Depleted starred: %d | Enriched starred: %d\n\n",
            sum(shap_df$depleted & shap_df$validated == "DOUBLY VALIDATED"),
            sum(shap_df$enriched)))

# Factor levels: top20 rank 1 at bottom (highest y), rank 20 at top (y=1)
# ggplot y-axis: levels[1] = y=1 (bottom), levels[20] = y=20 (top)
# We want rank 1 at TOP → levels should be reversed rank order
lev_labels <- rev(shap_df$feat_label)   # rank 20 first → rank 1 last (top of y-axis)
lev_feats  <- rev(shap_df$feature)

# y-axis label colours (parallel to lev_labels)
label_cols_vec <- label_color(lev_labels, lev_feats)

# ── Build long-format SHAP data for beeswarm ──────────────────────────────────
cat("Building beeswarm data...\n")
# Scale each feature's feature value (for color gradient)
X_scaled <- scale(X_shap[, top20_feat])  # 199 x 20

shap_long <- as.data.frame(shap_mat[, top20_feat]) %>%
  mutate(sample = seq_len(nrow(.))) %>%
  pivot_longer(-sample, names_to = "feature", values_to = "shap_value") %>%
  left_join(
    as.data.frame(X_scaled) %>%
      mutate(sample = seq_len(nrow(.))) %>%
      pivot_longer(-sample, names_to = "feature", values_to = "feat_val"),
    by = c("feature", "sample")
  ) %>%
  left_join(
    shap_df %>% select(feature, feat_label, validated, depleted),
    by = "feature"
  ) %>%
  mutate(
    feat_label = factor(feat_label, levels = lev_labels)
  )

set.seed(42)
shap_long$jitter_y <- as.numeric(shap_long$feat_label) +
                      runif(nrow(shap_long), -0.35, 0.35)

# ════════════════════════════════════════════════════════════════════════════════
cat("\n--- Figure 1: fig_shap_beeswarm.png ---\n")
# ════════════════════════════════════════════════════════════════════════════════
fig_beeswarm <- ggplot(shap_long,
    aes(x = shap_value, y = jitter_y, color = feat_val)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey30", linewidth = 0.5) +
  geom_point(alpha = 0.55, size = 1.2) +
  scale_color_gradient2(
    low = COL_DN, mid = "grey80", high = COL_UP, midpoint = 0,
    name = "Feature value\n(standardised)"
  ) +
  scale_y_continuous(breaks = seq_along(lev_labels), labels = lev_labels) +
  labs(
    title    = "SHAP beeswarm \u2014 eoCRC microbiome signature",
    subtitle = sprintf(
      "Combined RF | top-200 features | \u2605 = doubly validated (SHAP + MaAsLin2 v3) | %d features", n_val),
    x = "SHAP value  (positive \u2192 higher eoCRC probability)",
    y = NULL
  ) +
  theme_pres(11) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 9, color = label_cols_vec)
  )

ggsave("figures/fig_shap_beeswarm.png", fig_beeswarm, width = 11, height = 7, dpi = 300)
cat("  Saved: figures/fig_shap_beeswarm.png\n")
cat(sprintf("  Stars: %d | Red labels: %d | Blue labels: %d | Black labels: %d\n",
    sum(grepl("\u2605", lev_labels)),
    sum(label_cols_vec == COL_UP),
    sum(label_cols_vec == COL_DN),
    sum(label_cols_vec == COL_NS)))

# ════════════════════════════════════════════════════════════════════════════════
cat("\n--- Figure 2: fig_shap_bar.png ---\n")
# ════════════════════════════════════════════════════════════════════════════════
bar_df <- shap_df %>%
  mutate(
    feat_label_f = factor(feat_label, levels = lev_labels),
    fill_col = case_when(
      validated == "DOUBLY VALIDATED" & !depleted ~ "Enriched (validated)",
      validated == "DOUBLY VALIDATED" &  depleted ~ "Depleted (validated)",
      TRUE                                        ~ "SHAP only"
    )
  )

fill_pal <- c(
  "Enriched (validated)" = COL_UP,
  "Depleted (validated)" = COL_DN,
  "SHAP only"            = "grey60"
)

fig_bar <- ggplot(bar_df,
    aes(x = mean_abs_SHAP, y = feat_label_f, fill = fill_col)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = ifelse(validated == "DOUBLY VALIDATED", "\u2605", "")),
            hjust = -0.25, size = 4.5, color = "black") +
  scale_fill_manual(values = fill_pal, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title    = "Mean |SHAP| importance \u2014 top 20 features",
    subtitle = sprintf(
      "\u2605 = doubly validated (SHAP top-20 + MaAsLin2 v3 q<0.25) | %d features validated", n_val),
    x = "Mean |SHAP value|", y = NULL
  ) +
  theme_pres(11) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 9, color = label_cols_vec)
  )

ggsave("figures/fig_shap_bar.png", fig_bar, width = 9, height = 7, dpi = 300)
cat("  Saved: figures/fig_shap_bar.png\n")

# ════════════════════════════════════════════════════════════════════════════════
cat("\n--- Figure 3: fig_shap_dependence.png ---\n")
# ════════════════════════════════════════════════════════════════════════════════
top5_feat <- top20_feat[1:5]

dep_data <- lapply(seq_along(top5_feat), function(i) {
  fn  <- top5_feat[i]
  # Pick the feature with highest absolute Pearson correlation to this feature's SHAP
  # as the interaction color variable (mimic shapviz auto)
  other_feat <- setdiff(top20_feat, fn)
  cors <- sapply(other_feat, function(f2) abs(cor(shap_mat[, fn], X_shap[, f2])))
  int_feat <- other_feat[which.max(cors)]

  data.frame(
    feat_val  = X_shap[, fn],
    shap_val  = shap_mat[, fn],
    int_val   = X_shap[, int_feat],
    feature   = fn,
    int_name  = clean_sp(int_feat),
    clean_fn  = clean_sp(fn),
    rank      = i,
    stringsAsFactors = FALSE
  )
})

dep_plots <- lapply(dep_data, function(d) {
  lbl <- unique(d$clean_fn)
  int <- unique(d$int_name)
  validated_feat <- top5_feat[d$rank[1]] %in% V3_VALIDATED
  title_col <- if (validated_feat) {
    if (top5_feat[d$rank[1]] %in% V3_DEPLETED) COL_DN else COL_UP
  } else COL_NS

  ggplot(d, aes(x = feat_val, y = shap_val, color = int_val)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_point(alpha = 0.6, size = 1.4) +
    scale_color_gradient2(
      low = COL_DN, mid = "grey80", high = COL_UP, midpoint = 0,
      name = sprintf("Color:\n%s", int)
    ) +
    labs(
      title    = lbl,
      subtitle = if (validated_feat) "\u2605 doubly validated" else "SHAP only",
      x        = sprintf("%s (CLR-transformed)", lbl),
      y        = "SHAP value"
    ) +
    theme_pres(10) +
    theme(
      plot.title    = element_text(face = "bold.italic", size = 10, color = title_col),
      plot.subtitle = element_text(size = 8,
                                   color = if (validated_feat) title_col else "grey50"),
      legend.key.size = unit(0.35, "cm"),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 7)
    )
})

fig_dep <- wrap_plots(dep_plots, ncol = 2) +
  plot_annotation(
    title    = "SHAP dependence plots \u2014 top 5 features",
    subtitle = "SHAP value vs CLR-transformed abundance | color = most correlated interaction feature",
    theme    = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave("figures/fig_shap_dependence.png", fig_dep, width = 12, height = 10, dpi = 300)
cat("  Saved: figures/fig_shap_dependence.png\n")

# ════════════════════════════════════════════════════════════════════════════════
cat("\n--- Figure 4: fig_shap_combined.png (beeswarm + bar) ---\n")
# ════════════════════════════════════════════════════════════════════════════════
fig_combined <- (fig_beeswarm | fig_bar) +
  plot_layout(widths = c(1.2, 1)) +
  plot_annotation(
    title = "SHAP feature importance \u2014 eoCRC microbiome signature",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave("figures/fig_shap_combined.png", fig_combined, width = 18, height = 7, dpi = 300)
cat("  Saved: figures/fig_shap_combined.png\n")

# ════════════════════════════════════════════════════════════════════════════════
cat("\n--- Figure 5: fig_shap_panel.png (three-panel: beeswarm + bar + validation summary) ---\n")
# ════════════════════════════════════════════════════════════════════════════════
# Load LOSO fold AUROCs for the third panel (cross-model comparison)
folds <- readRDS("loso_primary_folds.rds")
folds_sp <- folds %>%
  filter(feature_set == "species") %>%
  select(study, RF_AUC, EN_AUC, XGB_AUC) %>%
  pivot_longer(-study, names_to = "model", values_to = "AUROC") %>%
  mutate(
    model = recode(model, RF_AUC = "RF", EN_AUC = "ElasticNet", XGB_AUC = "XGBoost"),
    model = factor(model, levels = c("RF", "XGBoost", "ElasticNet")),
    study_clean = gsub("_([0-9]{4}.*)", " \\1", study),
    study_clean = gsub("_", " ", study_clean)
  )

# Pooled means per model
pooled_means <- folds_sp %>%
  group_by(model) %>%
  summarise(mean_auc = mean(AUROC, na.rm = TRUE), .groups = "drop")

panel3 <- ggplot(folds_sp, aes(x = AUROC, y = reorder(study_clean, AUROC), color = model)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_point(size = 2.2, alpha = 0.8,
             position = position_dodge(width = 0.6)) +
  geom_vline(data = pooled_means,
             aes(xintercept = mean_auc, color = model),
             linetype = "dotted", linewidth = 1) +
  scale_color_manual(
    values = c("RF" = "#C0392B", "XGBoost" = "#2980B9", "ElasticNet" = "#27AE60"),
    name = "Model"
  ) +
  scale_x_continuous(limits = c(0.3, 1.0), breaks = seq(0.3, 1.0, 0.1)) +
  labs(
    title    = "Per-fold AUROC by model (species)",
    subtitle = "Dotted lines = pooled mean AUROC | Primary: eoCRC vs Young Controls",
    x = "AUROC", y = NULL
  ) +
  theme_pres(10) +
  theme(legend.position = "bottom")

# Compact beeswarm for panel (smaller)
fig_bee_compact <- fig_beeswarm +
  theme(legend.position = "none",
        plot.title = element_text(size = 11),
        plot.subtitle = element_text(size = 8))

fig_bar_compact <- fig_bar +
  theme(legend.position = "none",
        plot.title = element_text(size = 11),
        plot.subtitle = element_text(size = 8))

fig_panel <- (fig_bee_compact | fig_bar_compact | panel3) +
  plot_layout(widths = c(1.3, 1, 1)) +
  plot_annotation(
    title = "SHAP importance & cross-model AUROC \u2014 eoCRC microbiome signature",
    subtitle = sprintf(
      "Combined RF (top-200 features) | \u2605 = %d doubly validated (SHAP + MaAsLin2 v3)", n_val),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )

ggsave("figures/fig_shap_panel.png", fig_panel, width = 22, height = 8, dpi = 300)
cat("  Saved: figures/fig_shap_panel.png\n")

# ════════════════════════════════════════════════════════════════════════════════
cat("\n--- Figure 6: fig_shap_extended_combined.png (extended 200-feature beeswarm + bar) ---\n")
# ════════════════════════════════════════════════════════════════════════════════
# Beeswarm showing top-20 by mean |SHAP| (from 200-feature extended model)
# This is the main SHAP figure — same data as fig_beeswarm but labeled as "extended"

# Additional info in subtitle: note this is from the 200-feature combined RF
fig_bee_ext <- ggplot(shap_long,
    aes(x = shap_value, y = jitter_y, color = feat_val)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey30", linewidth = 0.5) +
  geom_point(alpha = 0.55, size = 1.2) +
  scale_color_gradient2(
    low = COL_DN, mid = "grey80", high = COL_UP, midpoint = 0,
    name = "Feature value\n(standardised)"
  ) +
  scale_y_continuous(breaks = seq_along(lev_labels), labels = lev_labels) +
  labs(
    title    = "SHAP beeswarm \u2014 eoCRC microbiome signature (extended model)",
    subtitle = sprintf(
      "Combined RF trained on top-200 features | Top 20 by mean |\u03a6| | \u2605 = %d doubly validated", n_val),
    x = "SHAP value  (positive \u2192 higher eoCRC probability)",
    y = NULL
  ) +
  theme_pres(11) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 9, color = label_cols_vec)
  )

fig_bar_ext <- ggplot(bar_df,
    aes(x = mean_abs_SHAP, y = feat_label_f, fill = fill_col)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = ifelse(validated == "DOUBLY VALIDATED", "\u2605", "")),
            hjust = -0.25, size = 4.5, color = "black") +
  scale_fill_manual(values = fill_pal, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title    = "Mean |\u03a6| importance \u2014 top 20 features (extended model)",
    subtitle = sprintf(
      "\u2605 = doubly validated | Red = enriched | Blue = depleted | %d features", n_val),
    x = "Mean |SHAP value|", y = NULL
  ) +
  theme_pres(11) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 9, color = label_cols_vec)
  )

fig_ext_combined <- (fig_bee_ext | fig_bar_ext) +
  plot_layout(widths = c(1.2, 1)) +
  plot_annotation(
    title = "SHAP feature importance \u2014 eoCRC microbiome signature (combined RF, 200 features)",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave("figures/fig_shap_extended_combined.png", fig_ext_combined,
       width = 18, height = 7, dpi = 300)
cat("  Saved: figures/fig_shap_extended_combined.png\n")

# ════════════════════════════════════════════════════════════════════════════════
cat("\n--- Figure 7: fig_shap_dependence_ext.png (dependence plots, top 5 extended) ---\n")
# ════════════════════════════════════════════════════════════════════════════════
# Same as fig_shap_dependence but explicitly from the extended (200-feature) model
# and using updated titles/annotations

dep_ext_plots <- lapply(dep_data, function(d) {
  lbl  <- unique(d$clean_fn)
  int  <- unique(d$int_name)
  feat <- top5_feat[d$rank[1]]
  validated_feat <- feat %in% V3_VALIDATED
  is_depleted    <- feat %in% V3_DEPLETED

  title_col <- if (validated_feat) {
    if (is_depleted) COL_DN else COL_UP
  } else COL_NS

  q_str <- if (validated_feat) {
    # Look up q-value from shap_df
    q <- shap_df$maaslin2_q[shap_df$feature == feat]
    if (!is.na(q)) sprintf("MaAsLin2 v3 q=%.4f", q) else ""
  } else ""

  dir_str <- if (validated_feat) {
    if (is_depleted) "\u2193 depleted in eoCRC" else "\u2191 enriched in eoCRC"
  } else "not in MaAsLin2 q<0.25"

  ggplot(d, aes(x = feat_val, y = shap_val, color = int_val)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_point(alpha = 0.6, size = 1.4) +
    scale_color_gradient2(
      low = COL_DN, mid = "grey80", high = COL_UP, midpoint = 0,
      name = sprintf("Color:\n%s", int)
    ) +
    labs(
      title    = if (validated_feat) paste0(lbl, "  \u2605") else lbl,
      subtitle = paste(dir_str, q_str, sep = " | "),
      x        = sprintf("%s (CLR-transformed)", lbl),
      y        = "SHAP value (\u03a6)"
    ) +
    theme_pres(10) +
    theme(
      plot.title    = element_text(face = "bold.italic", size = 10, color = title_col),
      plot.subtitle = element_text(size = 8,
                                   color = if (validated_feat) title_col else "grey50"),
      legend.key.size = unit(0.35, "cm"),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 7)
    )
})

fig_dep_ext <- wrap_plots(dep_ext_plots, ncol = 2) +
  plot_annotation(
    title    = "SHAP dependence plots \u2014 top 5 features (extended 200-feature model)",
    subtitle = sprintf(
      "SHAP value vs CLR abundance | color = interaction feature | \u2605 = doubly validated (%d / 5 shown)",
      sum(top5_feat %in% V3_VALIDATED)),
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )

ggsave("figures/fig_shap_dependence_ext.png", fig_dep_ext, width = 12, height = 10, dpi = 300)
cat("  Saved: figures/fig_shap_dependence_ext.png\n")

# ── Summary ────────────────────────────────────────────────────────────────────
cat("\n=== REGENERATION COMPLETE ===\n")
cat(sprintf("  v3 doubly validated in top-20: %d / 20\n", n_val))
cat(sprintf("  Enriched (red labels):   %d\n", sum(shap_df$enriched)))
cat(sprintf("  Depleted (blue labels):  %d\n", sum(shap_df$depleted & shap_df$validated == "DOUBLY VALIDATED")))
cat(sprintf("  Non-validated (black):   %d\n", 20 - n_val))
cat("\nNewly validated vs v2 (added 5 features):\n")
newly_v <- V3_VALIDATED[V3_VALIDATED %in% top20_feat]
old_val <- c("species:Peptostreptococcus stomatis","species:Gemella morbillorum",
             "species:Parvimonas micra","species:Dialister pneumosintes",
             "species:Alistipes finegoldii","species:Anaerotruncus colihominis",
             "species:Bacteroides caccae","species:[Clostridium] symbiosum")
new_adds <- setdiff(newly_v, old_val)
for (f in new_adds) cat(sprintf("  + %s\n", f))
cat("\nAll 7 figures saved to figures/ at 300 dpi.\n")
