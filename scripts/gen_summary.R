suppressPackageStartupMessages({
  library(dplyr)
})
setwd("/home/yugiu/eoCRC_analysis")

meta     <- readRDS("all_samples_metadata.rds") %>%
  mutate(grp = case_when(
    group=="CRC"     & age<50  ~ "eoCRC",
    group=="control" & age<50  ~ "young_ctrl",
    group=="CRC"     & age>=50 ~ "loCRC",
    group=="control" & age>=50 ~ "older_ctrl"))
pool_tbl <- readRDS("loso_pooled_auroc_summary.rds")
maas_all <- read.csv("maaslin2_results_v2/summaries/all_associations.csv",
                     stringsAsFactors=FALSE)
shap_df  <- read.csv("shap_top20_extended.csv", stringsAsFactors=FALSE)

pool_pri_sp <- pool_tbl %>% filter(analysis=="eoCRC_vs_YoungCtrl",
                                    feature_set=="species", model=="RF")
pool_pri_cb <- pool_tbl %>% filter(analysis=="eoCRC_vs_YoungCtrl",
                                    feature_set=="combined", model=="RF")
pool_sec_xg <- pool_tbl %>% filter(analysis=="eoCRC_vs_loCRC",
                                    feature_set=="species", model=="XGB")

n_eo <- sum(meta$grp=="eoCRC",      na.rm=TRUE)
n_yc <- sum(meta$grp=="young_ctrl", na.rm=TRUE)
n_lo <- sum(meta$grp=="loCRC",      na.rm=TRUE)
n_oc <- sum(meta$grp=="older_ctrl", na.rm=TRUE)

maas_eo_sp <- maas_all %>%
  filter(comparison=="eoCRC_vs_YoungCtrl", feature_set=="species", qval<0.25) %>%
  arrange(qval, desc(abs(coef)))
maas_lo_sp <- maas_all %>%
  filter(comparison=="loCRC_vs_OlderCtrl", feature_set=="species", qval<0.25) %>%
  arrange(qval, desc(abs(coef)))
maas_eo_pw <- maas_all %>%
  filter(comparison=="eoCRC_vs_YoungCtrl", feature_set=="pathways", qval<0.25) %>%
  arrange(qval)

mk_list <- function(df, label, n=10) {
  df10 <- head(df, n)
  sapply(seq_len(nrow(df10)), function(i) {
    dir <- if(df10$coef[i]>0) "up" else "down"
    nm  <- gsub("^species\\.","",df10$feature[i])
    nm  <- gsub("\\."," ",nm)
    nm  <- gsub("^[A-Z0-9_-]+  ","",nm)
    sprintf("  %d. %s (%s %s, q=%.4f)", i, nm, dir, label, df10$qval[i])
  })
}
top10_eo    <- mk_list(maas_eo_sp, "eoCRC")
top10_lo    <- mk_list(maas_lo_sp, "loCRC")
top10_eo_pa <- mk_list(maas_eo_pw, "eoCRC")

shap_rows <- sapply(seq_len(nrow(shap_df)), function(i) {
  q_str <- if(is.na(shap_df$maaslin2_q[i])) "NA" else sprintf("%.4f", shap_df$maaslin2_q[i])
  star  <- if(shap_df$validated[i]=="DOUBLY VALIDATED") "[star]" else ""
  nm    <- gsub("^species:","",shap_df$feature[i])
  sprintf("  %2d. %-45s | %s | q=%s %s",
          shap_df$rank[i], nm, shap_df$direction[i], q_str, star)
})

fn <- function(x, fmt) if(is.na(x)) "NA" else sprintf(fmt, x)

lines <- c(
  "# eoCRC Microbiome Meta-Analysis - Results Summary",
  "",
  paste0("**Date:** ", format(Sys.Date())),
  "**Institution:** University of Texas at Dallas (UTD)",
  "**Pipeline:** curatedMetagenomicData v3 -> MaAsLin2 -> LOSO-CV -> SHAP",
  "",
  "---",
  "",
  "## 1. Sample Counts",
  "",
  "| Group | N |",
  "|-------|---|",
  paste0("| eoCRC (CRC, age<50) | ", n_eo, " |"),
  paste0("| Young controls (healthy, age<50) | ", n_yc, " |"),
  paste0("| loCRC (CRC, age>=50) | ", n_lo, " |"),
  paste0("| Older controls (healthy, age>=50) | ", n_oc, " |"),
  paste0("| **Total** | **", n_eo+n_yc+n_lo+n_oc, "** |"),
  "",
  "**Feature matrices:** Species 1282x220 (post-filter); Pathways 1282x419 (unstratified)",
  "",
  "---",
  "",
  "## 2. LOSO-CV Classification Results",
  "",
  "### Primary: eoCRC vs Young Controls (8 valid folds, YuJ_2015 excluded)",
  "",
  "| Feature set | Model | Pooled AUROC | 95% CI | I2 | perm-p |",
  "|-------------|-------|-------------|--------|-----|--------|",
  paste0("| Species  | RF  | ",
         fn(pool_pri_sp$pooled_AUROC,"%.3f")," | ",
         fn(pool_pri_sp$CI_95_lo,"%.3f"),"-",fn(pool_pri_sp$CI_95_hi,"%.3f")," | ",
         fn(pool_pri_sp$I2_pct,"%.1f"),"% | ",fn(pool_pri_sp$perm_p,"%.4f")," |"),
  paste0("| Combined | RF  | ",
         fn(pool_pri_cb$pooled_AUROC,"%.3f")," | ",
         fn(pool_pri_cb$CI_95_lo,"%.3f"),"-",fn(pool_pri_cb$CI_95_hi,"%.3f")," | ",
         fn(pool_pri_cb$I2_pct,"%.1f"),"% | ",fn(pool_pri_cb$perm_p,"%.4f")," |"),
  "",
  "### Secondary: eoCRC vs loCRC (9 folds, loCRC downsampled in training)",
  "",
  "| Feature set | Model | Pooled AUROC | 95% CI | I2 | perm-p |",
  "|-------------|-------|-------------|--------|-----|--------|",
  paste0("| Species  | XGB | ",
         fn(pool_sec_xg$pooled_AUROC,"%.3f")," | ",
         fn(pool_sec_xg$CI_95_lo,"%.3f"),"-",fn(pool_sec_xg$CI_95_hi,"%.3f")," | ",
         fn(pool_sec_xg$I2_pct,"%.1f"),"% | ",fn(pool_sec_xg$perm_p,"%.4f")," |"),
  "",
  "---",
  "",
  "## 3. MaAsLin2 Differential Abundance",
  "",
  paste0("### eoCRC vs Young Controls - Species (q<0.25, total n=", nrow(maas_eo_sp), ")"),
  "",
  paste(top10_eo, collapse="\n"),
  "",
  paste0("### eoCRC vs Young Controls - Pathways (q<0.25, total n=", nrow(maas_eo_pw), ")"),
  "",
  paste(top10_eo_pa, collapse="\n"),
  "",
  paste0("### loCRC vs Older Controls - Species (q<0.25, total n=", nrow(maas_lo_sp), ")"),
  "",
  paste(top10_lo, collapse="\n"),
  "",
  "---",
  "",
  "## 4. SHAP Top 20 Features (Extended Combined RF, top-200)",
  "",
  "[star] = Doubly validated in both SHAP and MaAsLin2 (q<0.25)",
  "",
  paste(shap_rows, collapse="\n"),
  "",
  paste0("Doubly validated: ", sum(shap_df$validated=="DOUBLY VALIDATED"), "/20 features"),
  "",
  "---",
  "",
  "## 5. Key Biological Findings",
  "",
  "- eoCRC carries a detectable gut microbiome signature: species RF AUROC=0.704, perm-p=0.002",
  "- Signal is geographically consistent (I2=4.5%; 9 cohorts, Europe + East Asia)",
  "- Enriched in eoCRC: oral pathobionts Parvimonas micra, Gemella morbillorum, Peptostreptococcus stomatis, Dialister pneumosintes",
  "- Depleted in eoCRC: butyrate producers Intestinimonas butyriciproducens, Lachnospira eligens",
  "- Signal is primarily taxonomic (species RF p=0.002 vs pathway RF p=0.181 NS)",
  "- Oral-gut axis hypothesis: ectopic colonic colonization by oral pathobionts promotes tumorigenesis",
  "",
  "---",
  "",
  "## 6. Limitations",
  "",
  "1. Cross-sectional design - causality unestablished",
  "2. European/East Asian cohorts only - no Black/Hispanic representation",
  "3. Low eoCRC n=78 limits power for rare species",
  "4. No oral microbiome data - oral-gut axis inferred not measured",
  "5. Approximated SHAP (Monte Carlo nsim=30); direction may mislead for sparse features",
  "6. Pathway signal not significant (perm-p=0.181); functional interpretation limited",
  "",
  "---",
  "",
  "## 7. Future Directions",
  "",
  "1. Co-occurrence network analysis (SparCC) for polymicrobial biofilm structure",
  "2. Mendelian Randomization using FUT2/LCT genetic instruments",
  "3. Paired oral-stool study to directly test translocation of P. micra / G. morbillorum",
  "4. Dietary fiber intervention: test modulation of butyrate-producer:pathobiont ratio",
  "5. Validation in racially diverse US cohort (Parkland Health/UTSW)",
  "6. colibactin (pks island) quantification from raw metagenomic reads",
  "",
  "---",
  "",
  paste0("*Generated: ", format(Sys.time()), "*")
)

writeLines(lines, "summary.md")
cat("summary.md written:", file.size("summary.md"), "bytes\n")
