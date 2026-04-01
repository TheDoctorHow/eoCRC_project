#!/usr/bin/env Rscript
# Step 3b: Secondary permutation test + final pooled summary
# Fast version: pre-cache downsampled fold data (ntree=50, 4 cores)

suppressPackageStartupMessages({
  library(randomForest); library(pROC)
  library(compositions); library(meta)
  library(dplyr); library(parallel)
})
setwd("/home/yugiu/eoCRC_analysis")
set.seed(42)

N_PERMS    <- 1000
NTREE_PERM <- 50
N_CORES    <- max(1L, min(parallel::detectCores(logical=FALSE), 4L))
N_TOP      <- 50

cat(sprintf("Secondary permutation: %d perms, ntree=%d, %d cores\n\n",
            N_PERMS, NTREE_PERM, N_CORES))

# ── Load ───────────────────────────────────────────────────────────────────────
meta   <- readRDS("all_samples_metadata.rds") %>%
  mutate(maaslin_group = case_when(
    group == "CRC"     & age <  50 ~ "eoCRC",
    group == "CRC"     & age >= 50 ~ "loCRC",
    group == "control" & age <  50 ~ "young_ctrl",
    group == "control" & age >= 50 ~ "older_ctrl"
  ))
sp_mat <- readRDS("species_filtered.rds")
pa_mat <- readRDS("path_filtered.rds")

idx_s <- meta$maaslin_group %in% c("eoCRC", "loCRC")
sp_s  <- sp_mat[idx_s, ]; pa_s <- pa_mat[idx_s, ]
y_s   <- as.integer(meta$maaslin_group[idx_s] == "eoCRC")
st_s  <- meta$study_name[idx_s]

# ── Helpers (same as step3) ───────────────────────────────────────────────────
clr_mat <- function(X, pseudo=1e-6) as.matrix(compositions::clr(X + pseudo))

select_var_top <- function(X_tr, X_te, n=N_TOP) {
  top <- names(sort(apply(X_tr, 2, var), decreasing=TRUE))[seq_len(min(n, ncol(X_tr)))]
  list(tr=X_tr[, top, drop=FALSE], te=X_te[, top, drop=FALSE])
}

prep_fold <- function(sp_tr, sp_te, pa_tr, pa_te) {
  sc <- clr_mat(sp_tr); sct <- clr_mat(sp_te)
  pc <- clr_mat(pa_tr); pct <- clr_mat(pa_te)
  list(
    species  = select_var_top(sc,  sct),
    pathways = select_var_top(pc,  pct),
    combined = select_var_top(cbind(sc,pc), cbind(sct,pct))
  )
}

safe_auc <- function(truth, prob)
  tryCatch(as.numeric(pROC::auc(pROC::roc(truth, prob, quiet=TRUE))),
           error=function(e) NA_real_)

fit_rf_fast <- function(Xtr, ytr, Xte, yte, ntree=NTREE_PERM)
  tryCatch({
    m <- randomForest(x=Xtr, y=as.factor(ytr), ntree=ntree)
    safe_auc(yte, predict(m, Xte, type="prob")[,"1"])
  }, error=function(e) NA_real_)

# ── Build downsampled fold cache (same seeds as original run) ─────────────────
cat("Pre-computing downsampled fold cache...\n")
study_list <- unique(st_s)
cache <- list()
set.seed(42)    # match original run seed

for (s in study_list) {
  te_idx <- which(st_s == s)
  tr_idx <- which(st_s != s)
  y_te   <- y_s[te_idx]
  y_tr   <- y_s[tr_idx]

  # Downsample loCRC (class 0) to match eoCRC (class 1) n in training
  pos_tr  <- tr_idx[y_tr == 1]; neg_tr <- tr_idx[y_tr == 0]
  keep_neg <- sample(neg_tr, min(length(pos_tr), length(neg_tr)))
  dn_idx   <- c(pos_tr, keep_neg)
  y_dn     <- y_s[dn_idx]

  fsets <- prep_fold(sp_s[dn_idx,], sp_s[te_idx,],
                     pa_s[dn_idx,], pa_s[te_idx,])

  cache[[s]] <- list(dn_idx=dn_idx, te_idx=te_idx,
                     y_dn=y_dn, y_te=y_te, fsets=fsets)
  cat(sprintf("  %s: train=%d+%d  test=%d+%d\n",
              s, sum(y_dn==1), sum(y_dn==0),
              sum(y_te==1), sum(y_te==0)))
}
cat("\n")

# ── Permutation test ──────────────────────────────────────────────────────────
cat(sprintf("Running %d permutations on %d cores...\n", N_PERMS, N_CORES))

one_perm <- function(seed_i) {
  set.seed(seed_i)
  # Permute labels within each fold's downsampled training set
  out <- setNames(rep(NA_real_, 3), c("species","pathways","combined"))
  for (fs in names(out)) {
    fold_aucs <- vapply(names(cache), function(s) {
      fc    <- cache[[s]]
      y_te  <- fc$y_te
      if (length(unique(y_te)) < 2) return(NA_real_)
      y_perm <- sample(fc$y_dn)   # shuffle WITHIN downsampled training set
      fit_rf_fast(fc$fsets[[fs]]$tr, y_perm,
                  fc$fsets[[fs]]$te, y_te)
    }, numeric(1))
    out[fs] <- mean(fold_aucs, na.rm=TRUE)
  }
  out
}

perm_mat <- do.call(rbind,
  parallel::mclapply(seq_len(N_PERMS), one_perm,
                     mc.cores=N_CORES, mc.set.seed=FALSE))
secondary_perm <- as.data.frame(perm_mat)
saveRDS(secondary_perm, "loso_secondary_permutations.rds")
cat("Saved loso_secondary_permutations.rds\n\n")

# ── Reload primary results + build pooled summary table ───────────────────────
cat("Building pooled summary table...\n")

hanley_se <- function(auc, n1, n0) {
  q1 <- auc/(2-auc); q2 <- 2*auc^2/(1+auc)
  v  <- (auc*(1-auc) + (n1-1)*(q1-auc^2) + (n0-1)*(q2-auc^2)) / (n1*n0)
  sqrt(v) / (auc*(1-auc))
}

pool_dl <- function(aucs, n_cases, n_ctrls) {
  eps <- 1e-4
  ok  <- !is.na(aucs) & aucs > eps & aucs < 1-eps
  if (sum(ok) < 2)
    return(list(auc=mean(aucs,na.rm=TRUE), ci_lo=NA, ci_hi=NA,
                i2=NA, pval_het=NA, k=sum(ok)))
  a  <- aucs[ok]; n1 <- n_cases[ok]; n0 <- n_ctrls[ok]
  la <- log(a/(1-a))
  se <- mapply(hanley_se, a, n1, n0)
  res <- tryCatch(
    meta::metagen(TE=la, seTE=se, method.tau="DL", verbose=FALSE),
    error=function(e) NULL)
  if (is.null(res))
    return(list(auc=mean(a), ci_lo=NA, ci_hi=NA, i2=NA, pval_het=NA, k=sum(ok)))
  ilogit <- function(x) exp(x)/(1+exp(x))
  list(auc=ilogit(res$TE.random), ci_lo=ilogit(res$lower.random),
       ci_hi=ilogit(res$upper.random), i2=round(res$I2*100,1),
       pval_het=round(res$pval.Q,4), k=sum(ok))
}

make_table <- function(folds, perm_df, analysis) {
  rows <- list()
  for (fs in c("species","pathways","combined")) {
    df <- folds[folds$feature_set==fs,]
    for (mdl in c("RF","XGB","EN")) {
      aucs <- df[[paste0(mdl,"_AUC")]]
      pl   <- pool_dl(aucs, df$n_te_case, df$n_te_ctrl)
      obs  <- mean(aucs, na.rm=TRUE)
      p_p  <- if (mdl=="RF" && !is.null(perm_df))
                mean(perm_df[[fs]] >= obs, na.rm=TRUE) else NA_real_
      rows[[paste(fs,mdl)]] <- data.frame(
        analysis=analysis, feature_set=fs, model=mdl,
        k_folds=pl$k,
        mean_AUROC  =round(obs,      3),
        pooled_AUROC=round(pl$auc,   3),
        CI_95_lo    =round(pl$ci_lo, 3),
        CI_95_hi    =round(pl$ci_hi, 3),
        I2_pct      =pl$i2, p_heterog=pl$pval_het,
        perm_p      =if (!is.na(p_p)) round(p_p,4) else NA_real_,
        stringsAsFactors=FALSE)
    }
  }
  do.call(rbind, rows)
}

primary_folds   <- readRDS("loso_primary_folds.rds")
primary_perm    <- readRDS("loso_primary_permutations.rds")
secondary_folds <- readRDS("loso_secondary_folds.rds")

primary_pool   <- make_table(primary_folds,   primary_perm,   "eoCRC_vs_YoungCtrl")
secondary_pool <- make_table(secondary_folds, secondary_perm, "eoCRC_vs_loCRC")
all_pool       <- rbind(primary_pool, secondary_pool)

saveRDS(all_pool, "loso_pooled_auroc_summary.rds")
write.csv(all_pool, "loso_pooled_auroc_summary.csv", row.names=FALSE)

# ── Print results ─────────────────────────────────────────────────────────────
cat("\n\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║          POOLED AUROC SUMMARY — DerSimonian-Laird            ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("PRIMARY — eoCRC (n=78) vs Young Controls (n=121)  [8 valid folds]\n")
print(subset(all_pool, analysis=="eoCRC_vs_YoungCtrl", select=-analysis),
      row.names=FALSE)

cat("\nSECONDARY — eoCRC (n=78) vs loCRC (n=562)  [9 folds, loCRC downsampled]\n")
print(subset(all_pool, analysis=="eoCRC_vs_loCRC", select=-analysis),
      row.names=FALSE)

# Permutation detail
cat("\n\n── Permutation p-values (RF, 1000 perm) ──\n")
for (ana in c("eoCRC_vs_YoungCtrl","eoCRC_vs_loCRC")) {
  cat(sprintf("\n%s:\n", ana))
  perm_df  <- if (ana=="eoCRC_vs_YoungCtrl") primary_perm else secondary_perm
  folds_df <- if (ana=="eoCRC_vs_YoungCtrl") primary_folds else secondary_folds
  for (fs in c("species","pathways","combined")) {
    obs  <- mean(folds_df$RF_AUC[folds_df$feature_set==fs], na.rm=TRUE)
    pvec <- perm_df[[fs]]
    p    <- mean(pvec >= obs, na.rm=TRUE)
    q95  <- quantile(pvec, 0.95, na.rm=TRUE)
    cat(sprintf("  %-10s  obs=%.3f  95th-null=%.3f  perm-p=%.4f\n",
                fs, obs, q95, p))
  }
}

# Learning curve summary
lc <- readRDS("loso_primary_learning_curves.rds")
cat("\n\n── Learning curves — Primary (RF, mean AUROC across folds) ──\n")
lc_sum <- lc %>%
  group_by(feature_set, proportion) %>%
  summarise(mean_n=round(mean(n_train)),
            mean_AUROC=round(mean(RF_AUC,na.rm=TRUE),3),
            sd_AUROC  =round(sd(RF_AUC,  na.rm=TRUE),3),
            .groups="drop")
print(lc_sum, n=50)

# Fold-level detail (RF, species)
cat("\n\n── Per-fold RF AUROCs — Primary, species ──\n")
pf <- primary_folds[primary_folds$feature_set=="species",
                    c("study","n_te_case","n_te_ctrl","RF_AUC","EN_AUC")]
pf$RF_AUC <- round(pf$RF_AUC, 3); pf$EN_AUC <- round(pf$EN_AUC, 3)
print(pf, row.names=FALSE)

cat("\n── Per-fold RF AUROCs — Secondary, species (downsampled train) ──\n")
sf <- secondary_folds[secondary_folds$feature_set=="species",
                      c("study","n_te_case","n_te_ctrl","RF_AUC","EN_AUC")]
sf$RF_AUC <- round(sf$RF_AUC, 3); sf$EN_AUC <- round(sf$EN_AUC, 3)
print(sf, row.names=FALSE)

cat("\nStep 3 complete.\n")
