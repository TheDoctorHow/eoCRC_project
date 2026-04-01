#!/usr/bin/env Rscript
# Figure 4 — SHAP beeswarm (updated: 12 doubly validated features from v3 MaAsLin2)

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})

setwd("/home/yugiu/eoCRC_analysis")
dir.create("figures", showWarnings = FALSE)

# ── Palette & helpers (match step5_figures.R) ──────────────────────────────────
COL <- list(sig_up = "#C0392B", sig_dn = "#1A5276")

theme_pres <- function(base = 12)
  theme_bw(base_size = base) +
  theme(plot.title    = element_text(face = "bold", size = base + 1),
        plot.subtitle = element_text(size = base - 1, color = "grey40"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", color = NA),
        strip.text = element_text(face = "bold"))

clean_sp <- function(x) {
  x <- gsub("^species:", "", x)
  x <- gsub("[._]", " ", x)
  x <- gsub("sp  ", "sp. ", x)   # restore "sp." punctuation ("sp. X" encoded as "sp..X")
  ifelse(nchar(x) > 40, paste0(substr(x, 1, 38), "\u2026"), x)
}

# ── 13 doubly validated features (SHAP top-20 AND v3 MaAsLin2 q<0.25) ─────────
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
  "species:[Clostridium] symbiosum"   # q=0.2492 in v3
)

# Starred species depleted in eoCRC (MaAsLin2 v3 coef < 0): blue labels
V3_DEPLETED <- c(
  "species:Faecalibacterium prausnitzii",   # coef=-0.477
  "species:Lactobacillus rogosae",           # coef=-0.394
  "species:Turicimonas muris"               # coef=-0.513
)

# ── Load SHAP results ──────────────────────────────────────────────────────────
shap_ext   <- readRDS("shap_extended_results.rds")
shap_df    <- shap_ext$df
shap_mat   <- shap_ext$shap_mat
X_shap     <- shap_ext$Xsel
top20_feat <- shap_df$feature   # "species:X" format, length 20

cat(sprintf("Doubly validated in v3: %d / 20\n",
            sum(top20_feat %in% V3_VALIDATED)))
cat(sprintf("  Depleted starred: %d\n", sum(V3_VALIDATED %in% V3_DEPLETED)))

# ── Build long-format data for beeswarm ────────────────────────────────────────
shap_long <- as.data.frame(shap_mat) %>%
  mutate(sample = seq_len(nrow(.))) %>%
  pivot_longer(-sample, names_to = "feature", values_to = "shap_value") %>%
  inner_join(
    data.frame(
      feature  = colnames(X_shap),
      feat_val = as.numeric(scale(X_shap)[, 1:ncol(X_shap)]),
      stringsAsFactors = FALSE
    ) %>%
      group_by(feature) %>%
      mutate(rn = row_number()) %>%
      ungroup(),
    by = c("feature", "sample" = "rn")
  ) %>%
  filter(feature %in% top20_feat) %>%
  mutate(
    validated  = feature %in% V3_VALIDATED,
    depleted   = feature %in% V3_DEPLETED,
    clean_name = clean_sp(feature),
    feat_label = ifelse(validated,
                        paste0(clean_name, "  \u2605"),   # ★
                        clean_name),
    # Factor levels: rank order top→bottom on y-axis (top = highest y = rank 1)
    feat_label = factor(
      feat_label,
      levels = rev(ifelse(top20_feat %in% V3_VALIDATED,
                          paste0(clean_sp(top20_feat), "  \u2605"),
                          clean_sp(top20_feat)))
    )
  )

# ── Jitter y positions ─────────────────────────────────────────────────────────
set.seed(42)
shap_long$jitter_y <- as.numeric(shap_long$feat_label) +
                      runif(nrow(shap_long), -0.35, 0.35)

# ── y-axis label colours: starred enriched=red, starred depleted=blue, other=black
lev_labels <- levels(shap_long$feat_label)

# Map each label back to its feature to determine colour
label_to_feat <- setNames(shap_long$feature, shap_long$feat_label)
label_cols <- sapply(lev_labels, function(lbl) {
  feat <- label_to_feat[lbl]
  if (!grepl("\u2605", lbl)) return("black")          # non-validated: black
  if (feat %in% V3_DEPLETED)  return(COL$sig_dn)     # depleted starred: blue
  return(COL$sig_up)                                  # enriched starred: red
})

# ── Plot ───────────────────────────────────────────────────────────────────────
fig4 <- ggplot(shap_long, aes(x = shap_value, y = jitter_y, color = feat_val)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey30", linewidth = 0.5) +
  geom_point(alpha = 0.55, size = 1.2) +
  scale_color_gradient2(
    low = "#1A5276", mid = "grey80", high = "#C0392B",
    midpoint = 0, name = "Feature value\n(standardised)"
  ) +
  scale_y_continuous(
    breaks = seq_along(lev_labels),
    labels = lev_labels
  ) +
  labs(
    title    = "SHAP beeswarm — eoCRC microbiome signature",
    subtitle = "Combined RF | top-200 features | \u2605 = doubly validated (SHAP + MaAsLin2 v3) | 13 features",
    x = "SHAP value  (positive \u2192 higher eoCRC probability)",
    y = NULL
  ) +
  theme_pres(11) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 9, color = label_cols)
  )

# ── Save ───────────────────────────────────────────────────────────────────────
ggsave("figures/fig4_shap_beeswarm.png", fig4, width = 11, height = 7, dpi = 300)
ggsave("figures/fig4_shap_beeswarm.pdf", fig4, width = 11, height = 7)
cat("Saved: figures/fig4_shap_beeswarm.png + .pdf\n")
cat(sprintf("  Starred labels: %d\n", sum(grepl("\u2605", lev_labels))))
cat(sprintf("  Red (enriched): %d  |  Blue (depleted): %d  |  Black (non-validated): %d\n",
            sum(label_cols == COL$sig_up), sum(label_cols == COL$sig_dn),
            sum(label_cols == "black")))
cat("Done.\n")
