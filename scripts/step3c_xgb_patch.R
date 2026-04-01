#!/usr/bin/env Rscript
# Step 3c: Patch XGBoost AUROCs (xgboost v3 API fix) into existing fold results

suppressPackageStartupMessages({
  library(xgboost); library(pROC); library(compositions); library(meta); library(dplyr)
})
setwd("/home/yugiu/eoCRC_analysis")
set.seed(42)

N_TOP      <- 50
XGB_ROUNDS <- 100

meta <- readRDS("all_samples_metadata.rds") %>%
  mutate(maaslin_group = case_when(
    group == "CRC"     & age <  50 ~ "eoCRC",
    group == "CRC"     & age >= 50 ~ "loCRC",
    group == "control" & age <  50 ~ "young_ctrl",
    group == "control" & age >= 50 ~ "older_ctrl"
  ))
sp_mat <- readRDS("species_filtered.rds")
pa_mat <- readRDS("path_filtered.rds")

clr_mat <- function(X, pseudo=1e-6) as.matrix(compositions::clr(X + pseudo))
sel_top <- function(Xtr, Xte, n=N_TOP) {
  top <- names(sort(apply(Xtr,2,var), decreasing=TRUE))[seq_len(min(n,ncol(Xtr)))]
  list(tr=Xtr[,top,drop=FALSE], te=Xte[,top,drop=FALSE])
}
prep <- function(sp_tr, sp_te, pa_tr, pa_te) {
  sc  <- clr_mat(sp_tr); sct <- clr_mat(sp_te)
  pc  <- clr_mat(pa_tr); pct <- clr_mat(pa_te)
  list(species  = sel_top(sc,  sct),
       pathways = sel_top(pc,  pct),
       combined = sel_top(cbind(sc,pc), cbind(sct,pct)))
}
safe_auc <- function(truth, prob)
  tryCatch(as.numeric(pROC::auc(pROC::roc(truth, prob, quiet=TRUE))),
           error=function(e) NA_real_)
fit_xgb <- function(Xtr, ytr, Xte, yte) {
  tryCatch({
    dtr <- xgb.DMatrix(data=unname(Xtr), label=ytr)
    dte <- xgb.DMatrix(data=unname(Xte))
    m   <- xgb.train(params=list(objective="binary:logistic",
                                 eval_metric="auc", nthread=1),
                     data=dtr, nrounds=XGB_ROUNDS, verbose=0)
    safe_auc(yte, predict(m, dte))
  }, error=function(e) NA_real_)
}

run_xgb_loso <- function(sp, pa, y, studies, downsample=FALSE, label="") {
  rows <- list()
  for (s in unique(studies)) {
    te_idx <- which(studies==s); tr_idx <- which(studies!=s)
    y_te <- y[te_idx]; y_tr <- y[tr_idx]
    if (length(unique(y_te)) < 2) next
    if (downsample) {
      pos <- tr_idx[y_tr==1]; neg <- tr_idx[y_tr==0]
      tr_idx <- c(pos, sample(neg, min(length(pos),length(neg))))
      y_tr <- y[tr_idx]
    }
    fsets <- prep(sp[tr_idx,], sp[te_idx,], pa[tr_idx,], pa[te_idx,])
    for (fs in names(fsets)) {
      auc <- fit_xgb(fsets[[fs]]$tr, y_tr, fsets[[fs]]$te, y_te)
      cat(sprintf("  [%s|%-20s|%-9s] XGB=%.3f\n", label, s, fs, auc))
      rows[[paste(s,fs)]] <- data.frame(study=s, feature_set=fs,
                                         XGB_AUC=auc, stringsAsFactors=FALSE)
    }
  }
  do.call(rbind, rows)
}

# ── PRIMARY ──
cat("\n=== PRIMARY: eoCRC vs Young Controls ===\n")
idx_p <- meta$maaslin_group %in% c("eoCRC","young_ctrl")
xgb_p <- run_xgb_loso(sp_mat[idx_p,], pa_mat[idx_p,],
                       as.integer(meta$maaslin_group[idx_p]=="eoCRC"),
                       meta$study_name[idx_p], label="PRI")

pf <- readRDS("loso_primary_folds.rds")
pf$XGB_AUC <- NULL
pf <- merge(pf, xgb_p[,c("study","feature_set","XGB_AUC")],
            by=c("study","feature_set"), all.x=TRUE)
saveRDS(pf, "loso_primary_folds.rds")

# ── SECONDARY ──
cat("\n=== SECONDARY: eoCRC vs loCRC (downsampled) ===\n")
idx_s <- meta$maaslin_group %in% c("eoCRC","loCRC")
xgb_s <- run_xgb_loso(sp_mat[idx_s,], pa_mat[idx_s,],
                       as.integer(meta$maaslin_group[idx_s]=="eoCRC"),
                       meta$study_name[idx_s], downsample=TRUE, label="SEC")

sf <- readRDS("loso_secondary_folds.rds")
sf$XGB_AUC <- NULL
sf <- merge(sf, xgb_s[,c("study","feature_set","XGB_AUC")],
            by=c("study","feature_set"), all.x=TRUE)
saveRDS(sf, "loso_secondary_folds.rds")

# ── Rebuild final summary ──
cat("\n=== Rebuilding pooled summary ===\n")

hanley_se <- function(auc,n1,n0) {
  q1 <- auc/(2-auc); q2 <- 2*auc^2/(1+auc)
  v  <- (auc*(1-auc)+(n1-1)*(q1-auc^2)+(n0-1)*(q2-auc^2))/(n1*n0)
  sqrt(v)/(auc*(1-auc))
}
pool_dl <- function(aucs, n_cases, n_ctrls) {
  eps <- 1e-4; ok <- !is.na(aucs) & aucs>eps & aucs<1-eps
  if (sum(ok)<2) return(list(auc=mean(aucs,na.rm=TRUE),ci_lo=NA,ci_hi=NA,
                              i2=NA,pval_het=NA,k=sum(ok)))
  a <- aucs[ok]; n1 <- n_cases[ok]; n0 <- n_ctrls[ok]
  se <- mapply(hanley_se,a,n1,n0)
  res <- tryCatch(meta::metagen(TE=log(a/(1-a)), seTE=se,
                                method.tau="DL", verbose=FALSE),
                  error=function(e) NULL)
  if (is.null(res)) return(list(auc=mean(a),ci_lo=NA,ci_hi=NA,
                                 i2=NA,pval_het=NA,k=sum(ok)))
  il <- function(x) exp(x)/(1+exp(x))
  list(auc=il(res$TE.random), ci_lo=il(res$lower.random),
       ci_hi=il(res$upper.random), i2=round(res$I2*100,1),
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
        k_folds=pl$k, mean_AUROC=round(obs,3), pooled_AUROC=round(pl$auc,3),
        CI_95_lo=round(pl$ci_lo,3), CI_95_hi=round(pl$ci_hi,3),
        I2_pct=pl$i2, p_heterog=pl$pval_het,
        perm_p=if (!is.na(p_p)) round(p_p,4) else NA_real_,
        stringsAsFactors=FALSE)
    }
  }
  do.call(rbind, rows)
}

pf2   <- readRDS("loso_primary_folds.rds")
sf2   <- readRDS("loso_secondary_folds.rds")
pp    <- readRDS("loso_primary_permutations.rds")
sp_pm <- readRDS("loso_secondary_permutations.rds")

all_pool <- rbind(make_table(pf2, pp, "eoCRC_vs_YoungCtrl"),
                  make_table(sf2, sp_pm, "eoCRC_vs_loCRC"))
saveRDS(all_pool, "loso_pooled_auroc_summary.rds")
write.csv(all_pool, "loso_pooled_auroc_summary.csv", row.names=FALSE)

cat("\n\n╔═══════════════════════════════════════════════════════╗\n")
cat("║       FINAL POOLED AUROC SUMMARY (D-L, all models)    ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")
cat("PRIMARY — eoCRC vs Young Controls [8 folds]\n")
print(subset(all_pool, analysis=="eoCRC_vs_YoungCtrl", select=-analysis), row.names=FALSE)
cat("\nSECONDARY — eoCRC vs loCRC [9 folds, downsampled train]\n")
print(subset(all_pool, analysis=="eoCRC_vs_loCRC", select=-analysis), row.names=FALSE)
cat("\nStep 3c complete.\n")
