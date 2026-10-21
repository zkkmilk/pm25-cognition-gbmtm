library(nnet); library(mice); library(dplyr)
initial_options <- commandArgs(trailingOnly = FALSE); file_arg <- initial_options[grep("--file=", initial_options)]; if (length(file_arg) > 0) { SCRIPTS_DIR <- normalizePath(dirname(sub("--file=", "", file_arg[1]))) } else { SCRIPTS_DIR <- normalizePath(dirname(sys.frame(1)$ofile)) }
source(file.path(SCRIPTS_DIR, "config.R"))

log_msg("============================================================")
log_msg("Stage 3E [Sandbox]: Sensitivity (pp>=0.70 filter)")
log_msg("============================================================")

M <- STAGE3$mice_m; REF <- STAGE3$reference_group
anal_vars <- STAGE3$frozen_vars

bl <- read.csv(file.path(DIRS$data_processed, "baseline_with_region.csv"), stringsAsFactors = FALSE)
pp_mask <- bl$posterior_prob >= 0.70
n_full <- nrow(bl); n_filt <- sum(pp_mask)
log_msg(sprintf("pp>=0.70: %d/%d (%.1f%%)", n_filt, n_full, 100*n_filt/n_full))
for (g in sort(unique(bl$trajectory_group))) {
  log_msg(sprintf("  G%d: %d/%d (%.1f%%)", g,
    sum(bl$trajectory_group == g & pp_mask),
    sum(bl$trajectory_group == g),
    100*sum(bl$trajectory_group == g & pp_mask)/sum(bl$trajectory_group == g)))
}
bl_f <- bl[pp_mask, ]
stopifnot(length(unique(bl_f$trajectory_group)) >= 2)

log_msg("2. Missing mechanism (filtered)...")
miss_pct_f <- colMeans(is.na(bl_f[, anal_vars])) * 100
mice_vars_f <- names(miss_pct_f[miss_pct_f > 0])
log_msg(sprintf("  Variables with missing: %d", length(mice_vars_f)))

log_msg("3. MICE (filtered)...")
df_imp_f <- bl_f[, c(anal_vars, "trajectory_group")]
df_imp_f$trajectory_group <- factor(df_imp_f$trajectory_group)
if ("region" %in% names(df_imp_f)) {
  df_imp_f$region <- factor(df_imp_f$region, levels = c("东部","中部","西部","东北"))
}
meth_f <- rep("", ncol(df_imp_f)); names(meth_f) <- names(df_imp_f)
meth_f["trajectory_group"] <- ""
for (v in mice_vars_f) {
  if (v == "region") {
    meth_f[v] <- "polyreg"
  } else {
    uv <- unique(na.omit(df_imp_f[[v]]))
    meth_f[v] <- if (all(uv %in% c(0, 1))) "logreg" else "pmm"
  }
}
imp_f <- mice(df_imp_f, m = M, maxit = STAGE3$mice_maxit,
              method = meth_f, seed = STAGE3$mice_seed, print = FALSE)

imputed_f <- list()
for (m in 1:M) {
  ds <- complete(imp_f, m)
  ds$trajectory_group <- as.integer(as.character(ds$trajectory_group))
  if ("region" %in% names(ds)) ds$region <- as.character(ds$region)
  for (v in mice_vars_f) {
    if (v != "region") {
      uv <- unique(na.omit(bl_f[[v]]))
      if (all(uv %in% c(0, 1))) ds[[v]] <- round(pmax(pmin(ds[[v]], 1), 0))
    }
  }
  ds$trajectory_group <- factor(ds$trajectory_group)
  ds$trajectory_group <- relevel(ds$trajectory_group, ref = as.character(REF))
  imputed_f[[m]] <- ds
}
log_msg(sprintf("  %d datasets, %d rows", M, nrow(imputed_f[[1]])))

log_msg("4. VIF (filtered)...")
vif_vars_f <- setdiff(anal_vars, "region")
df_vif_f <- imputed_f[[1]][, vif_vars_f]
vif_f <- data.frame(variable = character(), VIF = numeric())
for (i in seq_along(vif_vars_f)) {
  v <- vif_vars_f[i]; preds <- setdiff(vif_vars_f, v)
  fit <- lm(as.formula(paste(v, "~", paste(preds, collapse = "+"))), data = df_vif_f)
  r2 <- summary(fit)$r.squared
  vif_f <- rbind(vif_f, data.frame(variable = v, VIF = round(1/(1-r2), 2)))
}
vif_f$variable <- sapply(vif_f$variable, label_var)
high_vif_f <- vif_f[vif_f$VIF > 5, ]
log_msg(sprintf("  VIF>5: %d", nrow(high_vif_f)))

log_msg("5. Regression (filtered)...")
mvars <- c("age","sex","education","marital_status","residence",
           "smoking","drinking","sleep_night","nap","employment","solid_fuel",
           "bmi","self_rated_health","adl_score","depression_cesd",
           "arthritis","cancer","liver_disease","kidney_disease","digestive_disease",
           "mobility_impairment","metabolic_burden","cardio_cerebrovascular","respiratory_burden")

BINARY_ALL_S <- c("sex","marital_status","residence","smoking","drinking",
  "employment","solid_fuel","arthritis","cancer","liver_disease","kidney_disease",
  "digestive_disease")
cont_vars_s <- setdiff(mvars, BINARY_ALL_S)
dummy_vars_s <- intersect(mvars, BINARY_ALL_S)

group_labels <- setdiff(levels(imputed_f[[1]]$trajectory_group), as.character(REF))

ds1 <- imputed_f[[1]]
X1c <- scale(ds1[, cont_vars_s, drop = FALSE])
fdat <- data.frame(X1c, ds1[, dummy_vars_s, drop = FALSE], trajectory_group = ds1$trajectory_group)
fit1 <- multinom(trajectory_group ~ ., data = fdat, trace = FALSE, maxit = 2000)
cm1 <- coef(fit1)
if (is.vector(cm1)) cm1 <- matrix(cm1, nrow = 1, dimnames = list(group_labels, names(cm1)))
coef_names <- colnames(cm1)

num_outcomes <- length(group_labels)
num_coefs <- length(coef_names)
n_tests_sens <- (length(coef_names) - 1) * num_outcomes
log_msg(sprintf("  Bonferroni family size (sensitivity) = %d", n_tests_sens))
coef_array <- array(NA, dim = c(num_outcomes, num_coefs, M),
                    dimnames = list(group_labels, coef_names, NULL))
se_array   <- array(NA, dim = c(num_outcomes, num_coefs, M),
                    dimnames = list(group_labels, coef_names, NULL))

for (m in 1:M) {
  ds <- imputed_f[[m]]
  X_cont_scaled <- scale(ds[, cont_vars_s, drop = FALSE])
  dat <- data.frame(X_cont_scaled,
                    ds[, dummy_vars_s, drop = FALSE],
                    trajectory_group = ds$trajectory_group)

  fit <- multinom(trajectory_group ~ ., data = dat, trace = FALSE, maxit = 2000)
  cm <- coef(fit)
  if (is.vector(cm)) cm <- matrix(cm, nrow = 1, dimnames = list(group_labels, names(cm)))
  se_cm <- summary(fit)$standard.errors
  if (is.vector(se_cm)) se_cm <- matrix(se_cm, nrow = 1, dimnames = list(group_labels, names(se_cm)))

  coef_array[,,m] <- cm[group_labels, coef_names]
  se_array[,,m]   <- se_cm[group_labels, coef_names]
}

pooled_coef <- apply(coef_array, c(1, 2), mean)
W_mat <- apply(se_array^2, c(1, 2), mean)
B_mat <- apply(coef_array, c(1, 2), var)
T_mat <- W_mat + (1 + 1/M) * B_mat
pooled_se <- sqrt(T_mat)

sens_rows <- list()
for (k in 1:num_outcomes) {
  for (cn in coef_names) {
    if (cn == "(Intercept)") next
    or_val <- exp(pooled_coef[k, cn])
    ci_lo  <- exp(pooled_coef[k, cn] - 1.96 * pooled_se[k, cn])
    ci_hi  <- exp(pooled_coef[k, cn] + 1.96 * pooled_se[k, cn])
    is_sep <- (pooled_se[k, cn] < 1e-6) || (or_val < 1e-4) || (ci_hi < 1e-4)
    if (is_sep) { p_val <- NA; p_adj <- NA; is_sig <- FALSE }
    else { z_val <- pooled_coef[k, cn] / (pooled_se[k, cn] + 1e-10)
           p_val <- 2 * (1 - pnorm(abs(z_val)))
           p_adj <- min(p_val * n_tests_sens, 1)
           is_sig <- (p_adj < 0.05) }

    label <- cn
    if (cn %in% names(VAR_LABELS)) {
      label <- VAR_LABELS[cn]
    }

    sens_rows[[length(sens_rows) + 1]] <- data.frame(
      outcome = group_labels[k], variable = label,
      OR = round(or_val, 3), CI_low = round(ci_lo, 3), CI_high = round(ci_hi, 3),
      p_value = if (is_sep) NA else round(p_val, 4),
      p_adj = if (is_sep) NA else p_adj,
      significant = if (is_sep) FALSE else is_sig, stringsAsFactors = FALSE)
  }
}
sens_df <- do.call(rbind, sens_rows)
log_msg(sprintf("  Significant (Bonferroni p_adj<0.05): %d/%d", sum(sens_df$significant), nrow(sens_df)))

log_msg("6. Comparing with main analysis...")
main_df <- read.csv(file.path(DIRS$tables, "multinomial_results.csv"))
main_m3 <- main_df[main_df$model == "Model3_全模型",
  c("outcome","variable","OR","CI_low","CI_high","p_value","p_adj")]
names(main_m3) <- c("outcome","variable","OR_main","CI_low_main","CI_high_main","p_main","p_adj_main")
merged <- merge(sens_df, main_m3, by = c("outcome","variable"))
merged$sig_main   <- !is.na(merged$p_adj_main) & (merged$p_adj_main < 0.05)
merged$direction_change <- (merged$OR > 1 & merged$OR_main < 1) | (merged$OR < 1 & merged$OR_main > 1)
merged$sig_lost   <- merged$sig_main & !merged$significant
merged$sig_gained <- !merged$sig_main & merged$significant
n_dir <- sum(merged$direction_change); n_lost <- sum(merged$sig_lost); n_gain <- sum(merged$sig_gained)
log_msg(sprintf("  Direction changes: %d, Sig lost: %d, Sig gained: %d (Bonferroni-matched)", n_dir, n_lost, n_gain))
write.csv(merged, file.path(DIRS$tables, "sensitivity_pp_filter.csv"), row.names = FALSE)

qa <- c("# Stage 3E QA: Sensitivity Analysis [Sandbox]",
        "", paste0("**Time**: ", format(Sys.time())),
        paste0("**Filter**: pp>=0.70, ", n_filt, "/", n_full),
        paste0("**Filtered VIF>5**: ", nrow(high_vif_f)),
        paste0("**Filtered M3 sig (Bonferroni)**: ", sum(sens_df$significant), "/", nrow(sens_df)),
        paste0("**Main M3 sig (Bonferroni)**: ", sum(merged$sig_main), "/", nrow(merged)),
        "", sprintf("Direction changes: %d", n_dir),
        sprintf("Sig lost: %d", n_lost),
        sprintf("Sig gained: %d", n_gain),
        "",
        "判定标准：与主分析同口径（Bonferroni，72 检验）。",
        "灵敏性子样本更小，部分主分析显著项跌出阈值属预期；关键判据为",
        "方向是否反转、以及是否出现主分析中不存在的新显著项。",
        "",
        if (n_dir == 0 && n_gain == 0)
          "### Conclusion: STABLE (no direction reversal, no new significance)"
        else "### Needs review")
writeLines(qa, file.path(DIRS$qc_reports, "stage3e_sensitivity_qc.md"))
log_msg("Stage 3E [Sandbox] complete")
