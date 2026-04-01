#!/usr/bin/env Rscript
# Figure 2 — LOSO-CV AUROC forest plot (clean rebuild from scratch)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

setwd("/home/yugiu/eoCRC_analysis")
dir.create("figures", showWarnings = FALSE)

# ── Fold data (species RF, eoCRC vs young controls; hardcoded from analysis) ──
# Listed top-to-bottom; y is assigned so top row = highest ggplot y value
fold_df <- data.frame(
  study = c("ThomasAM_2019_c", "VogtmannE_2016", "YachidaS_2019",
            "ZellerG_2014",    "HanniganGD_2017","FengQ_2015",
            "GuptaA_2019",     "WirbelJ_2018"),
  label = c("Thomas 2019 *",  "Vogtmann 2016", "Yachida 2019",
            "Zeller 2014",    "Hannigan 2017", "Feng 2015",
            "Gupta 2019",     "Wirbel 2018"),
  n1    = c(10,  8, 32, 3, 6, 4,  4,  7),   # test eoCRC cases
  n0    = c( 9, 11, 52, 8, 6, 3, 15, 17),   # test young-ctrl
  auc   = c(1.000, 0.932, 0.751, 0.750, 0.694, 0.667, 0.617, 0.504),
  stringsAsFactors = FALSE
) %>%
  mutate(
    # y increases upward in ggplot: Wirbel=1 (bottom), Thomas=8 (top fold row)
    y         = rev(seq_len(n())),
    row_label = sprintf("%s  (n=%d+%d)", label, n1, n0)
  )

# ── Hanley-McNeil 95% CI (logit-scale delta method) ──────────────────────────
hanley_ci <- function(auc, n1, n0) {
  if (auc >= 1 - 1e-6 || n1 < 2 || n0 < 2) return(c(NA_real_, NA_real_))
  auc <- max(1e-4, min(1 - 1e-4, auc))
  q1  <- auc / (2 - auc)
  q2  <- 2 * auc^2 / (1 + auc)
  se  <- sqrt((auc*(1-auc) + (n1-1)*(q1-auc^2) + (n0-1)*(q2-auc^2)) / (n1*n0))
  la  <- log(auc / (1 - auc))
  sel <- se / (auc * (1 - auc))
  il  <- function(x) exp(x) / (1 + exp(x))
  c(il(la - 1.96*sel), il(la + 1.96*sel))
}

fold_df <- fold_df %>%
  rowwise() %>%
  mutate(
    ci_lo = hanley_ci(auc, n1, n0)[1],
    ci_hi = hanley_ci(auc, n1, n0)[2],
    degen = auc >= 1 - 1e-6
  ) %>%
  ungroup()

# ── Pooled / meta-analysis summary ────────────────────────────────────────────
n_folds  <- nrow(fold_df)   # 8
pool_y   <- n_folds + 1.8   # 9.8 — pooled row sits above separator

POOL <- list(auc = 0.704, lo = 0.602, hi = 0.789, i2 = 4.5, p = 0.002)
SENS_AUC <- 0.704  # sensitivity analysis: excl. Thomas 2019

# Diamond vertices for pooled row
dmd <- data.frame(
  x = c(POOL$lo,  POOL$auc,      POOL$hi,  POOL$auc),
  y = c(pool_y,   pool_y + 0.45, pool_y,   pool_y - 0.45)
)

# ── Study colours ─────────────────────────────────────────────────────────────
COL_MAP <- c(
  ThomasAM_2019_c = "#2980B9", VogtmannE_2016  = "#C0392B",
  YachidaS_2019   = "#D35400", ZellerG_2014    = "#E91E63",
  HanniganGD_2017 = "#8E44AD", FengQ_2015      = "#E67E22",
  GuptaA_2019     = "#16A085", WirbelJ_2018    = "#27AE60"
)

# ── Layout constants ──────────────────────────────────────────────────────────
# x positions — all inside the panel (no clip="off" needed)
xlim_lo <- -0.36    # left panel edge; labels anchor at -0.34 so minimal gap
xlim_hi <-  1.50    # right panel edge; wide enough for "0.704 (0.602-0.789)"
x_label <- -0.34    # study name + n= label left anchor (hjust=0)
x_auc   <-  1.02    # AUROC value column left anchor (hjust=0)
x_annot <-  1.48    # I² / sensitivity right anchor (hjust=1)

# y positions
sep_y    <- n_folds + 1.0    # separator line y = 9.0
ref_y    <- 0.35             # "Chance" / "AUC=0.7" label y (below fold rows)
annot_y1 <- -0.50            # I² / perm-p annotation
annot_y2 <- -1.10            # sensitivity annotation
ylim_lo  <- -1.50
ylim_hi  <- pool_y + 0.9    # 10.7

# ── Build plot ────────────────────────────────────────────────────────────────
p <- ggplot(fold_df, aes(y = y)) +

  # Vertical reference lines (span full panel height)
  geom_vline(xintercept = 0.5, linetype = "dashed",
             color = "grey55", linewidth = 0.6) +
  geom_vline(xintercept = 0.7, linetype = "dotted",
             color = "#E74C3C", linewidth = 0.6) +

  # CI whiskers (non-degenerate folds only; Thomas 2019 omitted)
  geom_segment(
    data = filter(fold_df, !degen),
    aes(x = ci_lo, xend = ci_hi, yend = y),
    color = "grey35", linewidth = 0.85
  ) +

  # Fold dots (filled circles, coloured by study)
  geom_point(
    aes(x = auc, fill = study),
    shape = 21, size = 4.5, color = "white", show.legend = FALSE
  ) +
  scale_fill_manual(values = COL_MAP) +

  # Separator line between fold rows and pooled row
  geom_hline(yintercept = sep_y, color = "grey35", linewidth = 0.5) +

  # Pooled estimate diamond
  geom_polygon(
    data = dmd, aes(x = x, y = y),
    fill = "#2C3E50", color = NA, inherit.aes = FALSE
  ) +

  # ── Left column: study labels (inside panel) ──────────────────────────────
  geom_text(
    aes(label = row_label),
    x = x_label, hjust = 0, size = 3.6, color = "grey10"
  ) +
  annotate("text", x = x_label, y = pool_y, hjust = 0,
           size = 3.9, fontface = "bold", color = "grey10",
           label = "Pooled D-L  (k=8 folds)") +

  # ── Right column: AUROC values (inside panel) ─────────────────────────────
  geom_text(
    aes(label = sprintf("%.3f", auc)),
    x = x_auc, hjust = 0, size = 3.6, color = "grey10"
  ) +
  annotate("text", x = x_auc, y = pool_y, hjust = 0,
           size = 3.9, fontface = "bold", color = "grey10",
           label = sprintf("%.3f (%.3f-%.3f)", POOL$auc, POOL$lo, POOL$hi)) +

  # ── Reference line labels ─────────────────────────────────────────────────
  annotate("text", x = 0.5, y = ref_y, hjust = 0.5, size = 2.9,
           color = "grey45", label = "Chance") +
  annotate("text", x = 0.7, y = ref_y, hjust = 0.5, size = 2.9,
           color = "#E74C3C", label = "AUC=0.7") +

  # ── I² / perm-p annotation (bottom-right, inside panel) ──────────────────
  annotate("text", x = x_annot, y = annot_y1, hjust = 1, size = 3.1,
           color = "grey25",
           label = sprintf("I\u00B2=%.1f%%  |  perm-p=%.3f", POOL$i2, POOL$p)) +

  # ── Sensitivity annotation (inside panel) ────────────────────────────────
  annotate("text", x = x_annot, y = annot_y2, hjust = 1, size = 3.1,
           color = "grey25",
           label = sprintf("Sensitivity (excl. Thomas 2019): AUROC=%.3f", SENS_AUC)) +

  # ── Axes ─────────────────────────────────────────────────────────────────
  scale_x_continuous(
    breaks = seq(0.3, 1.0, by = 0.1),
    labels = sprintf("%.1f", seq(0.3, 1.0, by = 0.1))
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(
    xlim = c(xlim_lo, xlim_hi),
    ylim = c(ylim_lo, ylim_hi)
  ) +

  # ── Labels ───────────────────────────────────────────────────────────────
  labs(
    title    = "LOSO-CV AUROC forest plot: eoCRC vs Young Controls",
    subtitle = "Random Forest | Species features | Primary analysis",
    x        = "AUROC",
    y        = NULL,
    caption  = paste0(
      "* AUROC=1.0 with n=19 test samples; CI undefined for perfect classification; ",
      "D-L pooling down-weights this fold due to high variance"
    )
  ) +

  # ── Theme ─────────────────────────────────────────────────────────────────
  theme_bw(base_size = 12) +
  theme(
    plot.title         = element_text(face = "bold", size = 14),
    plot.subtitle      = element_text(size = 10, color = "grey40"),
    axis.title.x       = element_text(size = 12),
    axis.text.x        = element_text(size = 11),
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.caption       = element_text(size = 8.5, color = "grey40",
                                      face = "italic", hjust = 0,
                                      margin = margin(t = 4)),
    plot.margin        = margin(10, 120, 40, 10)
  )

# ── Save ─────────────────────────────────────────────────────────────────────
ggsave("figures/fig2_forest_plot.png", p, width = 10, height = 7, dpi = 300)
ggsave("figures/fig2_forest_plot.pdf", p, width = 10, height = 7)
cat("Saved: figures/fig2_forest_plot.png + .pdf\n")
cat("Done.\n")
