#!/usr/bin/env Rscript
# Figure 3 — Volcano plot (v3: proper CLR pre-computation)
# Reads from maaslin2_results_v3/summaries/all_associations.csv

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

setwd("/home/yugiu/eoCRC_analysis")
dir.create("figures", showWarnings = FALSE)

# ── Load v3 results ────────────────────────────────────────────────────────────
maas_v3 <- read.csv("maaslin2_results_v3/summaries/all_associations.csv",
                    stringsAsFactors = FALSE)

# ── Build volcano data for eoCRC vs young controls ─────────────────────────────
volc <- maas_v3 %>%
  filter(comparison == "eoCRC_vs_YoungCtrl", feature_set == "species") %>%
  mutate(
    neg_log_q = -log10(pmax(qval, 1e-10)),
    sig       = qval < 0.25,
    direction = case_when(
      sig & coef > 0 ~ "Enriched in eoCRC",
      sig & coef < 0 ~ "Depleted in eoCRC",
      TRUE           ~ "Not significant"
    ),
    # Clean species names
    label = feature,
    label = gsub("^species\\.", "", label),
    label = gsub("\\.", " ", label),
    # Fix double-space left by "sp. X" encoding (sp..X → "sp  X")
    label = gsub("sp  ", "sp. ", label),
    # Collapse any remaining multiple spaces and trim
    label = trimws(gsub("  +", " ", label))
  )

n_sig   <- sum(volc$sig, na.rm = TRUE)
n_up    <- sum(volc$sig & volc$coef > 0, na.rm = TRUE)
n_dn    <- sum(volc$sig & volc$coef < 0, na.rm = TRUE)

cat(sprintf("Significant (q<0.25): %d total  (%d UP, %d DOWN in eoCRC)\n",
            n_sig, n_up, n_dn))

# ── Species to always label ────────────────────────────────────────────────────
# Priority labels: biologically key + top by q-value
priority_label <- c(
  "Fusobacterium nucleatum",
  "Peptostreptococcus stomatis",
  "Parvimonas micra",
  "Gemella morbillorum",
  "Faecalibacterium prausnitzii",
  "Dialister pneumosintes",
  "Bacteroides cellulosilyticus",
  "Monoglobus pectinilyticus",
  "Anaerotruncus colihominis",
  "Romboutsia ilealis",
  "Bilophila wadsworthia",
  "Alistipes putredinis"
)

# Top 15 by q-value (significant only), then add priority hits
top_by_q <- volc %>%
  filter(sig) %>%
  arrange(qval, desc(abs(coef))) %>%
  head(15) %>%
  pull(label)

label_set <- union(top_by_q, priority_label)

KEY_TAXA <- c("Fusobacterium nucleatum", "Peptostreptococcus stomatis",
              "Parvimonas micra", "Gemella morbillorum",
              "Faecalibacterium prausnitzii")

LEG_LEVELS <- c("Key enriched taxa", "Key depleted taxa",
                 "Enriched in eoCRC", "Depleted in eoCRC", "Not significant")

volc <- volc %>%
  mutate(
    to_label = label %in% label_set & sig,
    key_taxa  = label %in% KEY_TAXA,
    legend_cat = factor(case_when(
      key_taxa & direction == "Enriched in eoCRC" ~ "Key enriched taxa",
      key_taxa & direction == "Depleted in eoCRC" ~ "Key depleted taxa",
      direction == "Enriched in eoCRC"            ~ "Enriched in eoCRC",
      direction == "Depleted in eoCRC"             ~ "Depleted in eoCRC",
      TRUE                                         ~ "Not significant"
    ), levels = LEG_LEVELS)
  )

# ── Colour / shape / size / alpha scales (5-level combined legend) ─────────────
LEG_COLOR <- c("Key enriched taxa" = "#C0392B", "Key depleted taxa" = "#1A5276",
               "Enriched in eoCRC" = "#C0392B", "Depleted in eoCRC" = "#1A5276",
               "Not significant"   = "#BDC3C7")
LEG_SHAPE <- c("Key enriched taxa" = 18, "Key depleted taxa" = 18,
               "Enriched in eoCRC" = 16, "Depleted in eoCRC" = 16,
               "Not significant"   = 16)
LEG_SIZE  <- c("Key enriched taxa" = 5.0, "Key depleted taxa" = 5.0,
               "Enriched in eoCRC" = 2.8, "Depleted in eoCRC" = 2.8,
               "Not significant"   = 1.4)
LEG_ALPHA <- c("Key enriched taxa" = 1.0,  "Key depleted taxa" = 1.0,
               "Enriched in eoCRC" = 0.72, "Depleted in eoCRC" = 0.72,
               "Not significant"   = 0.72)

# ── x-axis range: symmetric, clipped to meaningful range ───────────────────────
x_max <- ceiling(max(abs(volc$coef), na.rm = TRUE) * 10) / 10  # round up to 0.1
x_lim <- c(-x_max - 0.05, x_max + 0.05)

q_line_y  <- -log10(0.25)
annot_x   <- x_lim[2] - 0.02

# ── Plot ────────────────────────────────────────────────────────────────────────
fig3 <- ggplot(volc, aes(x = coef, y = neg_log_q,
                         color = legend_cat, shape = legend_cat)) +

  # q = 0.25 threshold line
  geom_hline(yintercept = q_line_y, linetype = "dashed",
             color = "grey50", linewidth = 0.7) +
  annotate("text", x = annot_x, y = q_line_y + 0.08,
           label = "q = 0.25", size = 3, color = "grey40", hjust = 1) +

  # All points — shape, size, alpha driven by legend_cat
  geom_point(aes(size = legend_cat, alpha = legend_cat)) +
  scale_color_manual(values = LEG_COLOR, name = NULL, breaks = LEG_LEVELS) +
  scale_shape_manual(values = LEG_SHAPE, name = NULL, breaks = LEG_LEVELS) +
  scale_size_manual( values = LEG_SIZE,  guide = "none") +
  scale_alpha_manual(values = LEG_ALPHA, guide = "none") +

  # Labels: most species (Fusobacterium handled separately to nudge it)
  geom_label_repel(
    data    = filter(volc, to_label & label != "Fusobacterium nucleatum"),
    aes(label = label),
    size               = 2.7,
    max.overlaps       = 40,
    box.padding        = 0.65,
    point.padding      = 0.35,
    label.padding      = 0.15,
    segment.color      = "grey40",
    segment.size       = 0.3,
    min.segment.length = 0.2,
    force              = 6,
    force_pull         = 0.5,
    seed               = 42,
    fill               = scales::alpha("white", 0.88),
    show.legend        = FALSE
  ) +

  # Fusobacterium nucleatum: nudged upward so it clears the q=0.25 line
  geom_label_repel(
    data    = filter(volc, label == "Fusobacterium nucleatum" & sig),
    aes(label = label),
    size               = 2.7,
    nudge_y            = 0.65,
    direction          = "y",
    box.padding        = 0.65,
    point.padding      = 0.35,
    label.padding      = 0.15,
    segment.color      = "grey40",
    segment.size       = 0.3,
    min.segment.length = 0.2,
    force              = 3,
    seed               = 42,
    fill               = scales::alpha("white", 0.88),
    show.legend        = FALSE
  ) +

  # Axes
  scale_x_continuous(
    breaks = seq(-2, 2, by = 0.5),
    labels = function(x) sprintf("%.1f", x)
  ) +
  coord_cartesian(xlim = x_lim) +

  # Labels
  labs(
    title    = "Differential abundance: eoCRC vs Young Controls",
    subtitle = sprintf(
      "MaAsLin2 | pre-computed CLR (compositions::clr, pseudocount=1e-6) | random effect: study | %d significant species (q<0.25)",
      n_sig),
    x = "Coefficient (CLR effect size)",
    y = expression(-log[10](q-value)),
    caption = paste0(
      "Coefficients are standardised (per SD of group indicator); multiply by ~2.05 for absolute CLR difference.\n",
      "Solid diamonds highlight key CRC-associated taxa. Positive coefficient = enriched in eoCRC."
    )
  ) +

  # Theme
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(size = 9, color = "grey35"),
    axis.title       = element_text(size = 12),
    axis.text        = element_text(size = 11),
    panel.grid.minor = element_blank(),
    legend.position  = c(0.15, 0.88),
    legend.background = element_rect(fill = scales::alpha("white", 0.85),
                                     color = "grey80", linewidth = 0.3),
    legend.text      = element_text(size = 10),
    legend.key.size  = unit(0.5, "lines"),
    plot.caption     = element_text(size = 8, color = "grey45", hjust = 0,
                                    margin = margin(t = 5))
  )

# ── Save ────────────────────────────────────────────────────────────────────────
ggsave("figures/fig3_volcano.png", fig3, width = 9, height = 7, dpi = 300)
ggsave("figures/fig3_volcano.pdf", fig3, width = 9, height = 7)
cat("Saved: figures/fig3_volcano.png + .pdf\n")
cat(sprintf("  UP eoCRC: %d  |  DOWN eoCRC: %d  |  Not significant: %d\n",
            n_up, n_dn, sum(!volc$sig, na.rm = TRUE)))
cat("Done.\n")
