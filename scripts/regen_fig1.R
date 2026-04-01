#!/usr/bin/env Rscript
# Regenerate Figure 1 (study design) — staggered x labels, caption inside plot

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})

setwd("/home/yugiu/eoCRC_analysis")
dir.create("figures", showWarnings = FALSE)

COL <- list(
  eoCRC      = "#E74C3C", young_ctrl = "#27AE60",
  loCRC      = "#2980B9", older_ctrl = "#8E44AD"
)

save_fig <- function(p, name, w=10, h=7) {
  ggsave(paste0("figures/",name,".png"), p, width=w, height=h, dpi=300)
  ggsave(paste0("figures/",name,".pdf"), p, width=w, height=h)
  cat(sprintf("  Saved: figures/%s.png + .pdf\n", name))
}

# ── Data ──────────────────────────────────────────────────────────────────────
meta <- readRDS("all_samples_metadata.rds") %>%
  mutate(grp = case_when(
    group=="CRC"     & age< 50 ~ "eoCRC",
    group=="CRC"     & age>=50 ~ "loCRC",
    group=="control" & age< 50 ~ "young_ctrl",
    group=="control" & age>=50 ~ "older_ctrl"
  ))

study_counts <- meta %>%
  filter(!is.na(grp)) %>%
  count(study_name, grp) %>%
  pivot_wider(names_from=grp, values_from=n, values_fill=0) %>%
  rename(study=study_name) %>%
  arrange(study) %>%
  mutate(total = eoCRC + young_ctrl + loCRC + older_ctrl,
         x     = row_number())

STUDIES <- study_counts$study
N_ST    <- nrow(study_counts)

short_name_map <- c(
  FengQ_2015       = "Feng 2015",
  GuptaA_2019      = "Gupta 2019",
  HanniganGD_2017  = "Hannigan 2017",
  ThomasAM_2019_c  = "Thomas 2019",
  VogtmannE_2016   = "Vogtmann 2016",
  WirbelJ_2018     = "Wirbel 2018",
  YachidaS_2019    = "Yachida 2019",
  YuJ_2015         = "Yu 2015",
  ZellerG_2014     = "Zeller 2014"
)
short_labels <- unname(short_name_map[STUDIES])

# ── Stacked rectangles (bar width = 0.6) ─────────────────────────────────────
seg_df <- study_counts %>%
  mutate(
    y_eo  = 0,          y_eo2 = eoCRC,
    y_yc  = eoCRC,      y_yc2 = eoCRC + young_ctrl,
    y_lo  = eoCRC + young_ctrl,
    y_lo2 = eoCRC + young_ctrl + loCRC,
    y_oc  = eoCRC + young_ctrl + loCRC,
    y_oc2 = total
  )

mk_rects <- function(df, y1col, y2col, fill, grp_lbl) {
  data.frame(xmin=df$x - 0.3, xmax=df$x + 0.3,
             ymin=df[[y1col]], ymax=df[[y2col]],
             fill=fill, grp=grp_lbl)
}
rect_df <- bind_rows(
  mk_rects(seg_df,"y_eo","y_eo2",  COL$eoCRC,      "eoCRC"),
  mk_rects(seg_df,"y_yc","y_yc2",  COL$young_ctrl, "Young ctrl"),
  mk_rects(seg_df,"y_lo","y_lo2",  COL$loCRC,      "loCRC"),
  mk_rects(seg_df,"y_oc","y_oc2",  COL$older_ctrl, "Older ctrl")
) %>% filter(ymin < ymax)

# ── Staggered x-axis labels (manually placed below y=0) ──────────────────────
# Odd-indexed cohorts: upper row (y = -16)
# Even-indexed cohorts: lower row (y = -34)
# Gap between bars (y=0) and upper row: 16 units — clearly separated
idx    <- seq_len(N_ST)
y_odd  <- -16
y_even <- -34

label_df <- data.frame(
  x     = idx,
  label = short_labels,
  y     = ifelse(idx %% 2 == 1, y_odd, y_even)
)

# Short tick segments from bar base (y=-2) down to just above each label row
tick_df <- data.frame(
  x    = idx,
  y0   = -2,
  y1   = ifelse(idx %% 2 == 1, y_odd + 4, y_even + 4)
)

# ── Key y positions ───────────────────────────────────────────────────────────
y_max_bar <- max(seg_df$y_oc2)
# Staggered labels sit at y_odd=-16, y_even=-34.
# LOSO arrow + label must clear the lower label row (y=-34) with room to spare.
loso_text <- -56     # LOSO label, below the lower stagger row
loso_line <- -68     # LOSO double-headed arrow, below the label

# ── Figure ────────────────────────────────────────────────────────────────────
fig1 <- ggplot() +

  # Stacked bars
  geom_rect(data=rect_df,
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=grp),
            color="white", linewidth=0.3) +
  scale_fill_manual(
    values = c("eoCRC"      = COL$eoCRC,
               "Young ctrl" = COL$young_ctrl,
               "loCRC"      = COL$loCRC,
               "Older ctrl" = COL$older_ctrl),
    name = "Group") +

  # Red dashed rectangle — example TEST fold (cohort 5)
  geom_rect(aes(xmin=4.65, xmax=5.35, ymin=-5, ymax=y_max_bar+14),
            fill=NA, color="#E74C3C", linetype="dashed", linewidth=1) +

  # Caption INSIDE plot near top-right of the dashed box
  annotate("text",
           x = 5.32, y = y_max_bar + 10,
           label = "example\ntest fold",
           size = 3.5, color = "#E74C3C", fontface = "italic",
           hjust = 1, vjust = 1) +

  # Tick segments connecting bars to staggered labels
  geom_segment(data=tick_df,
               aes(x=x, xend=x, y=y0, yend=y1),
               color="grey60", linewidth=0.35) +

  # Staggered x-axis labels
  geom_text(data=label_df,
            aes(x=x, y=y, label=label),
            size=3.8, hjust=0.5, vjust=1, color="grey15") +

  # LOSO label (above arrow)
  annotate("text",
           x=(N_ST+1)/2, y=loso_text,
           label="Leave-One-Study-Out Cross-Validation (LOSO-CV)",
           hjust=0.5, size=4.2, fontface="bold", color="grey20") +
  # LOSO double-headed arrow (below label)
  annotate("segment",
           x=0.7, xend=N_ST+0.3,
           y=loso_line, yend=loso_line,
           color="grey25", linewidth=0.9,
           arrow=arrow(ends="both", length=unit(5,"pt"), type="open")) +

  # Scales — suppress built-in x-axis text; y expanded for labels + LOSO
  scale_x_continuous(
    breaks = 1:N_ST,
    labels = NULL,
    expand = expansion(add = c(0.6, 0.6))
  ) +
  scale_y_continuous(
    expand = expansion(add = c(18, 30))  # 18 below loso_line (-68), 30 above bars
  ) +

  labs(
    title    = "Study design: 9-cohort LOSO-CV meta-analysis",
    subtitle = sprintf(
      "Total: %d samples  |  eoCRC=%d  |  Young ctrl=%d  |  loCRC=%d  |  Older ctrl=%d",
      sum(study_counts$total), sum(study_counts$eoCRC),
      sum(study_counts$young_ctrl), sum(study_counts$loCRC),
      sum(study_counts$older_ctrl)),
    x = NULL,
    y = "Sample count"
  ) +

  theme_bw(base_size=12) +
  theme(
    plot.title         = element_text(face="bold", size=16),
    plot.subtitle      = element_text(size=11, color="grey40"),
    axis.title.y       = element_text(size=14),
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    axis.text.y        = element_text(size=12),
    legend.text        = element_text(size=11),
    legend.title       = element_text(size=12, face="bold"),
    legend.position    = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.margin        = margin(t=8, r=10, b=30, l=8)  # generous bottom margin
  )

save_fig(fig1, "fig1_study_design", w=13, h=8)
cat("Done.\n")
