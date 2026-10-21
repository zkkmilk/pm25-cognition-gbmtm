library(nnet); library(dplyr); library(corrplot)
initial_options <- commandArgs(trailingOnly = FALSE); file_arg <- initial_options[grep("--file=", initial_options)]; if (length(file_arg) > 0) { SCRIPTS_DIR <- normalizePath(dirname(sub("--file=", "", file_arg[1]))) } else { SCRIPTS_DIR <- normalizePath(dirname(sys.frame(1)$ofile)) }
source(file.path(SCRIPTS_DIR, "config.R"))

log_msg("============================================================")
log_msg("Stage 3D: Post-imputation Diagnostics + Regression")
log_msg("============================================================")

M <- STAGE3$mice_m; REF <- STAGE3$reference_group
anal_vars <- STAGE3$frozen_vars

log_msg("1. Loading imputed datasets...")
imputed <- list()
for (m in 1:M) {
  ds <- read.csv(file.path(DIRS$imputed, sprintf("imputed_m%02d.csv", m)))
  ds$trajectory_group <- factor(ds$trajectory_group)
  ds$trajectory_group <- relevel(ds$trajectory_group, ref = as.character(REF))
  imputed[[m]] <- ds
}
log_msg(sprintf("Loaded %d datasets, %d rows each", M, nrow(imputed[[1]])))

log_msg("2. Post-imputation VIF...")
vif_vars <- setdiff(anal_vars, "region")
df_vif <- imputed[[1]][, vif_vars]
vif_result <- data.frame(variable = character(), VIF = numeric(), stringsAsFactors = FALSE)
for (i in seq_along(vif_vars)) {
  v <- vif_vars[i]
  preds <- setdiff(vif_vars, v)
  fmla <- as.formula(paste(v, "~", paste(preds, collapse = "+")))
  fit <- lm(fmla, data = df_vif)
  r2 <- summary(fit)$r.squared
  vif_result <- rbind(vif_result, data.frame(variable = v, VIF = round(1/(1-r2), 2)))
}
vif_result <- vif_result[order(-vif_result$VIF), ]
vif_result$variable <- sapply(vif_result$variable, label_var)
high_vif <- vif_result[vif_result$VIF > 5, ]
log_msg(sprintf("VIF>5: %d variables", nrow(high_vif)))
if (nrow(high_vif) > 0) for (i in 1:nrow(high_vif))
  log_msg(sprintf("  %s: VIF=%.2f", high_vif$variable[i], high_vif$VIF[i]))

log_msg("3. Correlation matrix...")
cor_mat <- cor(df_vif[, sapply(df_vif, is.numeric)], method = "pearson")
rownames(cor_mat) <- sapply(rownames(cor_mat), label_var)
colnames(cor_mat) <- sapply(colnames(cor_mat), label_var)
png(file.path(DIRS$figures, "stage3d_correlation_matrix.png"), width = 1600, height = 1400, res = 130)
corrplot(cor_mat, method = "color", type = "upper", tl.cex = 0.7,
         col = colorRampPalette(c("#D32F2F", "white", "#1976D2"))(200),
         addCoef.col = "black", number.cex = 0.5, diag = FALSE)
dev.off()

model1_vars <- c("age", "sex", "education", "marital_status", "residence")
model2_vars <- c(model1_vars, "smoking", "drinking", "sleep_night", "nap", "employment", "solid_fuel")
model3_vars <- c(model2_vars, "bmi", "self_rated_health", "adl_score", "depression_cesd",
  "arthritis", "cancer", "liver_disease", "kidney_disease", "digestive_disease",
  "mobility_impairment", "metabolic_burden", "cardio_cerebrovascular", "respiratory_burden")
models <- list("Model1_demographics" = model1_vars, "Model2_behavior_social" = model2_vars, "Model3_full model" = model3_vars)

BINARY_ALL <- c("sex","marital_status","residence","smoking","drinking",
  "employment","solid_fuel","arthritis","cancer","liver_disease","kidney_disease",
  "digestive_disease")

all_results <- list()
for (model_name in names(models)) {
  pred_vars <- models[[model_name]]
  log_msg(sprintf("\n--- %s (%d vars) ---", model_name, length(pred_vars)))

  cont_vars <- setdiff(pred_vars, BINARY_ALL)
  dummy_vars <- intersect(pred_vars, BINARY_ALL)

  group_labels <- setdiff(levels(imputed[[1]]$trajectory_group), as.character(REF))

  ds1 <- imputed[[1]]
  X1_cont <- scale(ds1[, cont_vars, drop = FALSE])
  fit1 <- multinom(trajectory_group ~ ., data = data.frame(
    X1_cont, ds1[, dummy_vars, drop = FALSE], trajectory_group = ds1$trajectory_group),
    trace = FALSE, maxit = 2000)
  cm1 <- coef(fit1)
  if (is.vector(cm1)) cm1 <- matrix(cm1, nrow = 1,
    dimnames = list(group_labels, names(cm1)))
  coef_names <- colnames(cm1)

  num_outcomes <- length(group_labels)
  num_coefs <- length(coef_names)
  coef_array <- array(NA, dim = c(num_outcomes, num_coefs, M),
                      dimnames = list(group_labels, coef_names, NULL))
  se_array   <- array(NA, dim = c(num_outcomes, num_coefs, M),
                      dimnames = list(group_labels, coef_names, NULL))

  for (m in 1:M) {
    ds <- imputed[[m]]
    X_cont_scaled <- scale(ds[, cont_vars, drop = FALSE])
    dat <- data.frame(X_cont_scaled,
                      ds[, dummy_vars, drop = FALSE],
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

  rows <- list()
  for (k in 1:num_outcomes) {
    for (cn in coef_names) {
      if (cn == "(Intercept)") next
      or_val <- exp(pooled_coef[k, cn])
      ci_lo  <- exp(pooled_coef[k, cn] - 1.96 * pooled_se[k, cn])
      ci_hi  <- exp(pooled_coef[k, cn] + 1.96 * pooled_se[k, cn])

      is_separated <- (pooled_se[k, cn] < 1e-6) || (or_val < 1e-4) || (ci_hi < 1e-4)

      if (is_separated) {
        p_val <- NA
        is_sig <- FALSE
      } else {
        z_val  <- pooled_coef[k, cn] / (pooled_se[k, cn] + 1e-10)
        p_val  <- 2 * (1 - pnorm(abs(z_val)))
        is_sig <- (p_val < 0.05)
      }

      label <- cn
      if (cn %in% names(VAR_LABELS)) {
        label <- VAR_LABELS[cn]
      }

      rows[[length(rows) + 1]] <- data.frame(
        model = model_name, outcome = group_labels[k], variable = label,
        OR = round(or_val, 3), CI_low = round(ci_lo, 3), CI_high = round(ci_hi, 3),
        p_value = if (is_separated) NA else round(p_val, 4),
        p_value_raw = if (is_separated) NA else p_val,
        significant = if (is_separated) FALSE else (p_val < 0.05), stringsAsFactors = FALSE)
    }
  }
  res_df <- do.call(rbind, rows)

  n_tests <- nrow(res_df)
  res_df$p_adj <- pmin(res_df$p_value_raw * n_tests, 1)
  res_df$significant <- !is.na(res_df$p_adj) & (res_df$p_adj < 0.05)

  all_results[[model_name]] <- res_df
  log_msg(sprintf("  Significant: %d/%d (uncorrected p<0.05) | %d/%d (Bonferroni p_adj<0.05, %d tests)",
                  sum(!is.na(res_df$p_value_raw) & res_df$p_value_raw < 0.05), n_tests,
                  sum(res_df$significant), n_tests, n_tests))
}

full_res <- do.call(rbind, all_results)
write.csv(full_res, file.path(DIRS$tables, "multinomial_results.csv"), row.names = FALSE)
write.csv(vif_result, file.path(DIRS$tables, "vif_post_imputation.csv"), row.names = FALSE)

log_msg("Generating forest plot with ggplot2...")
library(ggplot2)

forest_data <- all_results[["Model3_full model"]]

var_cat <- list(
  Demographics = c("Age","Sex","Education","Marital status","Residence","Employment"),
  Lifestyle    = c("Smoking","Drinking","Sleep night","Nap","Solid fuel use"),
  Physical     = c("BMI","Self-rated health","ADL score","Depression (CESD-10)","Mobility impairment"),
  Disease      = c("Arthritis","Cancer","Liver disease","Kidney disease","Digestive disease"),
  Composite    = c("Metabolic burden","Cardio-cerebrovascular","Respiratory burden"))

forest_data$cat <- NA
for (cn in names(var_cat)) {
  forest_data$cat[forest_data$variable %in% var_cat[[cn]]] <- cn
}

var_order <- rev(unlist(var_cat))
forest_data$variable <- factor(forest_data$variable, levels = var_order)

forest_data$plot_color <- ifelse(forest_data$significant, forest_data$cat, "Not significant (Bonferroni)")
forest_data$plot_color <- factor(forest_data$plot_color,
  levels = c(names(var_cat), "Not significant (Bonferroni)"))

cat_cols <- c(
  "Demographics"    = "#2196F3",
  "Lifestyle"       = "#4CAF50",
  "Physical"        = "#FF9800",
  "Disease"         = "#9C27B0",
  "Composite"       = "#F44336",
  "Not significant (Bonferroni)" = "gray65")

forest_data$outcome_label <- paste0("G", forest_data$outcome, " vs G", REF)

p_forest <- ggplot(forest_data, aes(x = OR, y = variable, color = plot_color)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.6) +
  geom_pointrange(aes(xmin = CI_low, xmax = CI_high),
    size = 0.5, linewidth = 0.8, shape = 19) +
  facet_wrap(~ outcome_label, nrow = 1) +
  scale_x_log10(
    breaks = c(0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0),
    labels = c("0.2", "0.5", "1.0", "2.0", "5.0", "10.0", "20.0")) +
  scale_y_discrete(expand = expansion(mult = c(0.02, 0.02))) +
  scale_color_manual(values = cat_cols, name = "Variable Categories") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(linetype = "dotted", color = "gray85"),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.8),
    strip.text = element_text(face = "bold", size = 14),
    axis.text.y = element_text(color = "black", size = 12),
    axis.text.x = element_text(color = "black", size = 12),
    axis.title.x = element_text(face = "bold", size = 13, margin = margin(t = 12)),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    legend.background = element_blank(),
    legend.box.margin = margin(t = -10),
    panel.border = element_rect(color = "black", linewidth = 0.8)) +
  labs(x = "Odds Ratio (95% CI) [Log Scale]", y = NULL) +
  guides(color = guide_legend(nrow = 1, title.position = "left", title.vjust = 0.5))

ggsave(file.path(DIRS$figures, "stage3d_forest_plot.png"),
  plot = p_forest, width = 14, height = 15, dpi = 300)

qa <- c("# Stage 3D QA: Post-imputation Diagnostics + Regression",
        "", paste0("**Time**: ", format(Sys.time())),
        "", "## VIF (post-imputation)",
        sprintf("VIF > 5: %d variables", nrow(high_vif)),
        "", "## Model Summary",
        sprintf("- M1: %d/%d significant", sum(all_results[["Model1_demographics"]]$significant), nrow(all_results[["Model1_demographics"]])),
        sprintf("- M2: %d/%d significant", sum(all_results[["Model2_behavior_social"]]$significant), nrow(all_results[["Model2_behavior_social"]])),
        sprintf("- M3: %d/%d significant", sum(all_results[["Model3_full model"]]$significant), nrow(all_results[["Model3_full model"]])),
        "", paste0("![Forest](../figures/stage3d_forest_plot.png)"),
        "", "> Stage 3D complete.")
writeLines(qa, file.path(DIRS$qc_reports, "stage3d_regression_qc.md"))

log_msg("Stage 3D complete")
