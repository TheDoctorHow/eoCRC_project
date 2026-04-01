#!/usr/bin/env Rscript
# Step 3: LOSO-CV ML Classification
#
# PRIMARY:   eoCRC (78) vs young controls (121) — 8 valid folds (YuJ_2015 skipped)
# SECONDARY: eoCRC (78) vs loCRC (562)          — 9 folds, loCRC downsampled in train
#
# Feature sets: species | pathways | combined (each CLR'd within fold, top-50 by var)
# Models: Random Forest (ntree=500) | XGBoost (nrounds=100) | ElasticNet (alpha=0.5)
# Pooling: DerSimonian-Laird on logit-AUROC with Hanley-McNeil SE
# Permutation: 1000 shuffles, RF only (ntree=100), parallel

suppressPackageStartupMessages({
  library(randomForest)
  library(xgboost)
  library(glmnet)
  library(pROC)
  library(meta)
  library(compositions)
  library(dplyr)
  library(parallel)
})

setwd("/home/yugiu/eoCRC_analysis")
set.seed(42)

# ── Parameters ────────────────────────────────────────────────────────────────
N_TOP      <- 50       # features after variance selection
NTREE      <- 500
NTREE_PERM <- 100      # faster trees for permutation null
XGB_ROUNDS <- 100
N_PERMS    <- 1000
LC_PROPS   <- c(0.25, 0.50, 0.75, 1.00)
N_CORES    <- max(1L, min(parallel::detectCores(logical = FALSE), 4L))

cat(sprintf("Cores available for permutation: %d\n\n", N_CORES))

# ── 1. Load data ──────────────────────────────────────────────────────────────
cat("Loading data...\n")
meta   <- readRDS("all_samples_metadata.rds")
sp_mat <- readRDS("species_filtered.rds")   # 1282 x 220, values 0-100
pa_mat <- readRDS("path_filtered.rds")      # 1282 x 419

meta <- meta %>%
  mutate(maaslin_group = case_when(
    group == "CRC"     & age <  50 ~ "eoCRC",
    group == "CRC"     & age >= 50 ~ "loCRC",
    group == "control" & age <  50 ~ "young_ctrl",
    group == "control" & age >= 50 ~ "older_ctrl"
  ))

stopifnot(all(rownames(sp_mat) == meta$sample_id))
stopifnot(all(rownames(pa_mat) == meta$sample_id))

# ── 2. Helper functions ───────────────────────────────────────────────────────

# CLR applied independently to each matrix (row-wise; scale-invariant so no /100 needed)
clr_mat <- function(X, pseudo = 1e-6) {
  as.matrix(compositions::clr(X + pseudo))
}

# Variance-based feature selection: learn cutoff from training only
select_var_top <- function(X_tr, X_te, n = N_TOP) {
  vars <- apply(X_tr, 2, var)
  top  <- names(sort(vars, decreasing = TRUE))[seq_len(min(n, ncol(X_tr)))]
  list(tr = X_tr[, top, drop = FALSE],
       te = X_te[, top, drop = FALSE])
}

# Prepare all 3 feature sets for one fold (CLR within fold, then select)
prep_fold <- function(sp_tr, sp_te, pa_tr, pa_te) {
  sp_c <- clr_mat(sp_tr);  sp_ct <- clr_mat(sp_te)
  pa_c <- clr_mat(pa_tr);  pa_ct <- clr_mat(pa_te)
  cb_tr <- cbind(sp_c, pa_c);   cb_te <- cbind(sp_ct, pa_ct)

  list(
    species  = select_var_top(sp_c,  sp_ct),
    pathways = select_var_top(pa_c,  pa_ct),
    combined = select_var_top(cb_tr, cb_te)
  )
}

# ── Model wrappers ────────────────────────────────────────────────────────────
safe_auc <- function(truth, prob) {
  tryCatch(
    as.numeric(pROC::auc(pROC::roc(truth, prob, quiet = TRUE))),
    error = function(e) NA_real_
  )
}

fit_rf <- function(Xtr, ytr, Xte, yte, ntree = NTREE) {
  tryCatch({
    m <- randomForest(x = Xtr, y = as.factor(ytr), ntree = ntree)
    safe_auc(yte, predict(m, Xte, type = "prob")[, "1"])
  }, error = function(e) NA_real_)
}

fit_xgb <- function(Xtr, ytr, Xte, yte) {
  tryCatch({
    # xgboost v3 API: use xgb.DMatrix + xgb.train; strip col names first
    Xtr2 <- unname(Xtr); Xte2 <- unname(Xte)
    dtr  <- xgboost::xgb.DMatrix(data = Xtr2, label = ytr)
    dte  <- xgboost::xgb.DMatrix(data = Xte2)
    m    <- xgboost::xgb.train(
               params  = list(objective = "binary:logistic",
                              eval_metric = "auc", nthread = 1),
               data    = dtr, nrounds = XGB_ROUNDS, verbose = 0)
    safe_auc(yte, predict(m, dte))
  }, error = function(e) NA_real_)
}

fit_en <- function(Xtr, ytr, Xte, yte) {
  tryCatch({
    nf <- max(3L, min(5L, as.integer(min(table(ytr)))))
    m  <- glmnet::cv.glmnet(Xtr, ytr, alpha = 0.5,
                             family = "binomial", nfolds = nf)
    p  <- as.numeric(predict(m, Xte, s = "lambda.min", type = "response"))
    safe_auc(yte, p)
  }, error = function(e) NA_real_)
}

# ── 3. LOSO-CV ─────────────────────────────────────────────────────────────────
loso_cv <- function(sp, pa, y, studies,
                    label     = "",
                    downsample = FALSE,
                    verbose    = TRUE) {
  study_list <- unique(studies)
  rows <- list()

  for (s in study_list) {
    te_idx <- which(studies == s)
    tr_idx <- which(studies != s)
    y_te   <- y[te_idx]
    y_tr   <- y[tr_idx]

    if (length(unique(y_te)) < 2) {
      if (verbose) cat(sprintf("  SKIP fold '%s': single class in test set\n", s))
      next
    }

    # Downsample majority class in TRAINING only
    if (downsample) {
      pos_tr <- tr_idx[y_tr == 1]; neg_tr <- tr_idx[y_tr == 0]
      keep_neg <- sample(neg_tr, min(length(pos_tr), length(neg_tr)))
      tr_idx   <- c(pos_tr, keep_neg)
      y_tr     <- y[tr_idx]
    }

    fsets <- prep_fold(sp[tr_idx, ], sp[te_idx, ],
                       pa[tr_idx, ], pa[te_idx, ])

    for (fs in names(fsets)) {
      Xtr <- fsets[[fs]]$tr;  Xte <- fsets[[fs]]$te

      auc_rf  <- fit_rf(Xtr, y_tr, Xte, y_te)
      auc_xgb <- fit_xgb(Xtr, y_tr, Xte, y_te)
      auc_en  <- fit_en(Xtr, y_tr, Xte, y_te)

      if (verbose) cat(sprintf(
        "  [%s|%-10s|%-9s] RF=%.3f XGB=%.3f EN=%.3f  tr=%d+%d  te=%d+%d\n",
        label, s, fs,
        auc_rf, auc_xgb, auc_en,
        sum(y_tr == 1), sum(y_tr == 0),
        sum(y_te == 1), sum(y_te == 0)
      ))

      rows[[paste(s, fs, sep = "__")]] <- data.frame(
        study      = s, feature_set = fs,
        n_tr_case  = sum(y_tr == 1), n_tr_ctrl = sum(y_tr == 0),
        n_te_case  = sum(y_te == 1), n_te_ctrl = sum(y_te == 0),
        RF_AUC     = auc_rf,
        XGB_AUC    = auc_xgb,
        EN_AUC     = auc_en,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# ── 4. DerSimonian-Laird pooling ──────────────────────────────────────────────
# SE via Hanley-McNeil formula on AUROC scale, then delta-method to logit scale
hanley_mcneil_se <- function(auc, n1, n0) {
  q1  <- auc / (2 - auc)
  q2  <- 2 * auc^2 / (1 + auc)
  v   <- (auc * (1 - auc) +
          (n1 - 1) * (q1 - auc^2) +
          (n0 - 1) * (q2 - auc^2)) / (n1 * n0)
  # delta method: SE(logit(AUC)) = SE(AUC) / (AUC*(1-AUC))
  sqrt(v) / (auc * (1 - auc))
}

pool_dl <- function(aucs, n_cases, n_ctrls) {
  eps  <- 1e-4
  ok   <- !is.na(aucs) & aucs > eps & aucs < 1 - eps
  if (sum(ok) < 2)
    return(list(auc = mean(aucs, na.rm=TRUE), ci_lo=NA, ci_hi=NA,
                i2=NA, pval_het=NA, k=sum(ok)))

  a <- aucs[ok]; n1 <- n_cases[ok]; n0 <- n_ctrls[ok]
  logit_a <- log(a / (1 - a))
  se      <- mapply(hanley_mcneil_se, a, n1, n0)

  res <- tryCatch(
    meta::metagen(TE = logit_a, seTE = se, method.tau = "DL", verbose = FALSE),
    error = function(e) NULL
  )
  if (is.null(res))
    return(list(auc = mean(a), ci_lo=NA, ci_hi=NA, i2=NA, pval_het=NA, k=sum(ok)))

  ilogit <- function(x) exp(x) / (1 + exp(x))
  list(
    auc      = ilogit(res$TE.random),
    ci_lo    = ilogit(res$lower.random),
    ci_hi    = ilogit(res$upper.random),
    i2       = round(res$I2 * 100, 1),
    pval_het = round(res$pval.Q, 4),
    k        = sum(ok)
  )
}

# ── 5. Permutation test (RF, ntree=NTREE_PERM) ───────────────────────────────
perm_test <- function(sp, pa, y, studies, n_perm = N_PERMS, label = "") {
  study_list <- unique(studies)

  # Pre-compute CLR + feature selection once (both are y-independent)
  cat(sprintf("  Building fold cache for permutation test (%s)...\n", label))
  cache        <- list()
  valid_studies <- c()

  for (s in study_list) {
    te_idx <- which(studies == s)
    tr_idx <- which(studies != s)
    if (length(unique(y[te_idx])) < 2) next
    valid_studies <- c(valid_studies, s)
    cache[[s]] <- list(
      tr_idx = tr_idx, te_idx = te_idx,
      fsets  = prep_fold(sp[tr_idx, ], sp[te_idx, ],
                         pa[tr_idx, ], pa[te_idx, ])
    )
  }
  cat(sprintf("  Valid folds: %d | %d perms on %d cores\n",
              length(valid_studies), n_perm, N_CORES))

  one_perm <- function(seed_i) {
    set.seed(seed_i)
    y_p <- sample(y)   # global shuffle preserving study structure

    out <- setNames(rep(NA_real_, 3), c("species", "pathways", "combined"))
    for (fs in names(out)) {
      fold_aucs <- vapply(valid_studies, function(s) {
        fc   <- cache[[s]]
        y_tr <- y_p[fc$tr_idx]; y_te <- y_p[fc$te_idx]
        if (length(unique(y_te)) < 2) return(NA_real_)
        fit_rf(fc$fsets[[fs]]$tr, y_tr,
               fc$fsets[[fs]]$te, y_te, ntree = NTREE_PERM)
      }, numeric(1))
      out[fs] <- mean(fold_aucs, na.rm = TRUE)
    }
    out
  }

  perm_mat <- do.call(rbind,
    parallel::mclapply(seq_len(n_perm), one_perm, mc.cores = N_CORES,
                       mc.set.seed = FALSE)
  )
  as.data.frame(perm_mat)
}

# ── 6. Learning curves (RF, stratified subsampling) ───────────────────────────
learning_curves <- function(sp, pa, y, studies, props = LC_PROPS, label = "") {
  study_list <- unique(studies)
  rows <- list()

  for (s in study_list) {
    te_idx <- which(studies == s)
    tr_idx <- which(studies != s)
    y_te   <- y[te_idx]
    if (length(unique(y_te)) < 2) next

    pos_all <- tr_idx[y[tr_idx] == 1]
    neg_all <- tr_idx[y[tr_idx] == 0]

    for (prop in props) {
      n_pos <- max(2L, round(length(pos_all) * prop))
      n_neg <- max(2L, round(length(neg_all) * prop))
      sub   <- c(sample(pos_all, min(n_pos, length(pos_all))),
                 sample(neg_all, min(n_neg, length(neg_all))))
      y_tr  <- y[sub]
      if (length(unique(y_tr)) < 2) next

      fsets <- prep_fold(sp[sub, ], sp[te_idx, ],
                         pa[sub, ], pa[te_idx, ])

      for (fs in names(fsets)) {
        auc <- fit_rf(fsets[[fs]]$tr, y_tr, fsets[[fs]]$te, y_te)
        rows[[paste(s, prop, fs)]] <- data.frame(
          study      = s, proportion = prop,
          n_train    = length(sub),
          n_tr_case  = sum(y_tr == 1), n_tr_ctrl = sum(y_tr == 0),
          feature_set = fs, RF_AUC = auc,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

# ── 7. Pooling summary helper ─────────────────────────────────────────────────
make_pool_table <- function(folds, perm_df, analysis) {
  rows <- list()
  for (fs in c("species", "pathways", "combined")) {
    df <- folds[folds$feature_set == fs, ]
    for (mdl in c("RF", "XGB", "EN")) {
      aucs <- df[[paste0(mdl, "_AUC")]]
      pl   <- pool_dl(aucs, df$n_te_case, df$n_te_ctrl)
      obs  <- mean(aucs, na.rm = TRUE)

      p_perm <- if (mdl == "RF" && !is.null(perm_df))
        mean(perm_df[[fs]] >= obs, na.rm = TRUE) else NA_real_

      rows[[paste(fs, mdl)]] <- data.frame(
        analysis    = analysis, feature_set = fs, model = mdl,
        k_folds     = pl$k,
        mean_AUROC  = round(obs,       3),
        pooled_AUROC= round(pl$auc,    3),
        CI_95_lo    = round(pl$ci_lo,  3),
        CI_95_hi    = round(pl$ci_hi,  3),
        I2_pct      = pl$i2,
        p_heterog   = pl$pval_het,
        perm_p      = if (!is.na(p_perm)) round(p_perm, 4) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# ══════════════════════════════════════════════════════════════════════════════
# PRIMARY ANALYSIS: eoCRC vs Young Controls
# ══════════════════════════════════════════════════════════════════════════════
cat("\n\n══ PRIMARY: eoCRC vs Young Controls ══\n\n")

idx_p  <- meta$maaslin_group %in% c("eoCRC", "young_ctrl")
sp_p   <- sp_mat[idx_p, ];  pa_p   <- pa_mat[idx_p, ]
y_p    <- as.integer(meta$maaslin_group[idx_p] == "eoCRC")
st_p   <- meta$study_name[idx_p]
cat(sprintf("n=%d  eoCRC=%d  young_ctrl=%d\n\n", sum(idx_p), sum(y_p), sum(!y_p)))

cat("--- LOSO-CV ---\n")
primary_folds <- loso_cv(sp_p, pa_p, y_p, st_p, label = "PRI")
saveRDS(primary_folds, "loso_primary_folds.rds")

cat("\n--- Learning curves (RF) ---\n")
primary_lc <- learning_curves(sp_p, pa_p, y_p, st_p, label = "PRI")
saveRDS(primary_lc, "loso_primary_learning_curves.rds")
cat(sprintf("  Saved %d curve rows\n", nrow(primary_lc)))

cat("\n--- Permutation test ---\n")
primary_perm <- perm_test(sp_p, pa_p, y_p, st_p, n_perm = N_PERMS, label = "PRI")
saveRDS(primary_perm, "loso_primary_permutations.rds")

# ══════════════════════════════════════════════════════════════════════════════
# SECONDARY ANALYSIS: eoCRC vs loCRC (downsampled training)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n\n══ SECONDARY: eoCRC vs loCRC (loCRC downsampled in train) ══\n\n")

idx_s  <- meta$maaslin_group %in% c("eoCRC", "loCRC")
sp_s   <- sp_mat[idx_s, ];  pa_s   <- pa_mat[idx_s, ]
y_s    <- as.integer(meta$maaslin_group[idx_s] == "eoCRC")
st_s   <- meta$study_name[idx_s]
cat(sprintf("n=%d  eoCRC=%d  loCRC=%d\n\n", sum(idx_s), sum(y_s), sum(!y_s)))

cat("--- LOSO-CV ---\n")
secondary_folds <- loso_cv(sp_s, pa_s, y_s, st_s,
                            label = "SEC", downsample = TRUE)
saveRDS(secondary_folds, "loso_secondary_folds.rds")

cat("\n--- Permutation test ---\n")
secondary_perm <- perm_test(sp_s, pa_s, y_s, st_s,
                             n_perm = N_PERMS, label = "SEC")
saveRDS(secondary_perm, "loso_secondary_permutations.rds")

# ══════════════════════════════════════════════════════════════════════════════
# POOL + SUMMARISE
# ══════════════════════════════════════════════════════════════════════════════
cat("\n\n══ POOLED AUROC SUMMARY ══\n\n")

primary_pool   <- make_pool_table(primary_folds,   primary_perm,   "eoCRC_vs_YoungCtrl")
secondary_pool <- make_pool_table(secondary_folds, secondary_perm, "eoCRC_vs_loCRC")
all_pool       <- rbind(primary_pool, secondary_pool)

saveRDS(all_pool, "loso_pooled_auroc_summary.rds")
write.csv(all_pool, "loso_pooled_auroc_summary.csv", row.names = FALSE)

# Pretty-print
cat("PRIMARY — eoCRC vs Young Controls\n")
print(subset(all_pool, analysis == "eoCRC_vs_YoungCtrl",
             select = -analysis), row.names = FALSE)

cat("\nSECONDARY — eoCRC vs loCRC\n")
print(subset(all_pool, analysis == "eoCRC_vs_loCRC",
             select = -analysis), row.names = FALSE)

# Permutation test detail (RF)
cat("\n\n── Empirical permutation p-values (RF only) ──\n")
for (ana in c("eoCRC_vs_YoungCtrl", "eoCRC_vs_loCRC")) {
  cat(sprintf("\n%s:\n", ana))
  perm_df <- if (ana == "eoCRC_vs_YoungCtrl") primary_perm else secondary_perm
  folds_df <- if (ana == "eoCRC_vs_YoungCtrl") primary_folds else secondary_folds
  for (fs in c("species", "pathways", "combined")) {
    obs  <- mean(folds_df$RF_AUC[folds_df$feature_set == fs], na.rm = TRUE)
    pvec <- perm_df[[fs]]
    p    <- mean(pvec >= obs, na.rm = TRUE)
    q    <- quantile(pvec, 0.95, na.rm = TRUE)
    cat(sprintf("  %-10s  obs=%.3f  95th-perm=%.3f  perm-p=%.4f\n",
                fs, obs, q, p))
  }
}

# Learning curve summary
cat("\n\n── Learning curves (RF mean AUROC across folds) ──\n")
lc_sum <- primary_lc %>%
  group_by(feature_set, proportion) %>%
  summarise(
    mean_n_train = round(mean(n_train)),
    mean_AUROC   = round(mean(RF_AUC, na.rm=TRUE), 3),
    sd_AUROC     = round(sd(RF_AUC,   na.rm=TRUE), 3),
    .groups = "drop"
  )
print(lc_sum, n = 50)

cat("\nStep 3 complete.\n")
