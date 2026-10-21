library(gbmt)
library(doParallel)
library(dplyr)

initial_options <- commandArgs(trailingOnly = FALSE); file_arg <- initial_options[grep("--file=", initial_options)]; if (length(file_arg) > 0) { SCRIPTS_DIR <- normalizePath(dirname(sub("--file=", "", file_arg[1]))) } else { SCRIPTS_DIR <- normalizePath(dirname(sys.frame(1)$ofile)) }
source(file.path(SCRIPTS_DIR, "config.R"))

calc_entropy <- function(posterior) {
  ng <- ncol(posterior)
  n  <- nrow(posterior)
  eps <- 1e-16
  p <- pmax(pmin(posterior, 1 - eps), eps)
  1 - sum(-rowSums(p * log(p))) / (n * log(ng))
}

log_msg("============================================================")
log_msg("Stage 2: GBMTM joint trajectory modeling")
log_msg("============================================================")

log_msg("1. loading trajectory data...")
sim_data <- read.csv(file.path(DIRS$data_processed, "trajectory_long.csv"))
sim_data$ID   <- as.character(sim_data$ID)
sim_data$year <- as.numeric(sim_data$year)
log_msg(sprintf("   individuals: %d, rows: %d", n_distinct(sim_data$ID), nrow(sim_data)))

log_msg(sprintf("2. parallel GBMTM fitting (d=%d, scaling=%d, ng=%s, cores=%d)...",
  STAGE2$d, STAGE2$scaling, paste(STAGE2$ng_range, collapse=","), STAGE2$n_cores))

cl <- makeCluster(STAGE2$n_cores)
registerDoParallel(cl)
clusterExport(cl, c("sim_data", "STAGE2", "DIRS"))

model_results <- foreach(ng = STAGE2$ng_range, .packages = "gbmt") %dopar% {
  start_t <- Sys.time()
  model <- gbmt(
    x.names = c("cog_score", "PM25"),
    unit    = "ID",
    time    = "time",
    d       = STAGE2$d,
    ng      = ng,
    data    = sim_data,
    scaling = STAGE2$scaling,
    pruning = STAGE2$pruning
  )
  run_time <- round(as.numeric(difftime(Sys.time(), start_t, units = "secs")), 2)
  rds_path <- file.path(DIRS$models, sprintf("gbmtm_ng%d.rds", ng))
  saveRDS(model, rds_path)
  list(ng = ng, run_time = run_time,
       aic = round(AIC(model), 2), bic = round(BIC(model), 2), model = model)
}
stopCluster(cl)
log_msg("   all model fits completed")

log_msg("3. computing fit indices...")
compare_list <- list()
for (res in model_results) {
  ng <- res$ng; mod <- res$model
  post <- mod$posterior
  ek <- calc_entropy(post)
  appa_m <- round(mean(mod$appa), 4)
  gp <- prop.table(table(mod$assign))
  min_prop <- round(min(gp), 4)
  conv <- mod$em[nrow(mod$em), "converged"] == 1
  compare_list[[paste0("ng", ng)]] <- data.frame(
    ng = ng, AIC = res$aic, BIC = res$bic,
    entropy = ek, APPA_mean = appa_m,
    min_group_prop = min_prop,
    converged = as.integer(conv),
    iterations = nrow(mod$em),
    run_time_sec = res$run_time,
    stringsAsFactors = FALSE
  )
  log_msg(sprintf("   ng=%d AIC=%.1f BIC=%.1f Ek=%.4f APPA=%.4f conv=%d",
    ng, res$aic, res$bic, ek, appa_m, as.integer(conv)))
}
compare_table <- do.call(rbind, compare_list)

log_msg("4. model selection...")
conv_models <- compare_table[compare_table$converged == 1, ]
stopifnot(nrow(conv_models) > 0)
best_idx <- which.min(conv_models$BIC)
best_ng  <- conv_models$ng[best_idx]
log_msg(sprintf("   lowest BIC BIC: ng=%d", best_ng))

sel <- STAGE2$selection
if (conv_models$entropy[best_idx] < sel$min_entropy ||
    conv_models$APPA_mean[best_idx] < sel$min_appa ||
    conv_models$min_group_prop[best_idx] < sel$min_group_prop) {
  for (i in order(conv_models$ng, decreasing = TRUE)) {
    r <- conv_models[i, ]
    if (r$entropy >= sel$min_entropy &&
        r$APPA_mean >= sel$min_appa &&
        r$min_group_prop >= sel$min_group_prop) {
      best_ng <- r$ng
      log_msg(sprintf("   >> falling back to ng=%d", best_ng))
      break
    }
  }
}
log_msg(sprintf("   Final: ng=%d", best_ng))
best_model <- model_results[[which(sapply(model_results, `[[`, "ng") == best_ng)]]$model

write.csv(compare_table, file.path(DIRS$tables, "model_comparison.csv"), row.names = FALSE)

log_msg("6. group relabeling...")
ng <- best_ng
fitted_vals <- best_model$fitted
group_stats <- data.frame(group = 1:ng, mean_pm25 = NA, cog_slope = NA, mean_cog = NA)
for (g in 1:ng) {
  fit <- fitted_vals[[g]]
  group_stats$mean_pm25[g] <- mean(fit[, 2], na.rm = TRUE)
  group_stats$cog_slope[g]  <- fit[nrow(fit), 1] - fit[1, 1]
  group_stats$mean_cog[g]   <- mean(fit[, 1], na.rm = TRUE)
}
group_stats$risk_score <- scale(group_stats$mean_pm25) - scale(group_stats$cog_slope)
group_stats$new_label  <- rank(group_stats$risk_score, ties.method = "first")

for (i in 1:nrow(group_stats)) {
  log_msg(sprintf("   old=%d -> new=G%d  PM25=%.1f  slope=%.3f",
    group_stats$group[i], group_stats$new_label[i],
    group_stats$mean_pm25[i], group_stats$cog_slope[i]))
}

log_msg("7. saving trajectory assignment...")
old_labels  <- best_model$assign
new_labels  <- setNames(group_stats$new_label, group_stats$group)
mapped_assign <- new_labels[as.character(old_labels)]

post_df <- as.data.frame(best_model$posterior)
old_for_new <- match(1:ng, group_stats$new_label)
post_df <- post_df[, old_for_new]
names(post_df) <- paste0("prob_G", 1:ng)

assign_df <- data.frame(
  ID = names(old_labels),
  trajectory_group = as.integer(mapped_assign),
  posterior_prob   = post_df[cbind(1:nrow(post_df), as.integer(mapped_assign))],
  stringsAsFactors = FALSE
)
assign_df <- cbind(assign_df, post_df)

write.csv(assign_df, file.path(DIRS$models, "final_trajectory_assignment.csv"), row.names = FALSE)
write.csv(group_stats, file.path(DIRS$models, "group_reordering.csv"), row.names = FALSE)
log_msg(sprintf("   individuals: %d", nrow(assign_df)))

log_msg("8. generating trajectory plot...")
time_pts <- c(0, 2, 4, 7)
time_labels <- c("0(2011)", "2(2013)", "4(2015)", "7(2018)")
cols <- c("#2196F3", "#4CAF50", "#FF9800", "#F44336")[1:ng]

png(file.path(DIRS$figures, "stage2_trajectory_by_group.png"), width = 1200, height = 560, res = 120)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

yl <- range(sapply(fitted_vals, function(x) x[, 1]))
plot(time_pts, fitted_vals[[1]][, 1], type = "n", ylim = yl,
     xlab = "Time", ylab = "Cognition (raw)", main = "Cognition", xaxt = "n")
axis(1, at = time_pts, labels = time_labels)
abline(h = 0, lty = 3, col = "gray")
for (g in 1:ng) {
  fit <- fitted_vals[[old_for_new[g]]]
  lines(time_pts, fit[, 1], col = cols[g], lwd = 2.5, type = "o", pch = 19)
}
legend("bottomleft", paste0("G", 1:ng), col = cols[1:ng], lwd = 2, cex = 0.8)

yl <- range(sapply(fitted_vals, function(x) x[, 2]))
plot(time_pts, fitted_vals[[1]][, 2], type = "n", ylim = yl,
     xlab = "Time", ylab = "PM2.5 (raw)", main = "PM2.5", xaxt = "n")
axis(1, at = time_pts, labels = time_labels)
abline(h = 0, lty = 3, col = "gray")
for (g in 1:ng) {
  fit <- fitted_vals[[old_for_new[g]]]
  lines(time_pts, fit[, 2], col = cols[g], lwd = 2.5, type = "o", pch = 19)
}
legend("topright", paste0("G", 1:ng), col = cols[1:ng], lwd = 2, cex = 0.8)
dev.off()

png(file.path(DIRS$figures, "stage2_posterior_histogram.png"), width = 1000, height = 800, res = 130)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
for (g in 1:ng) {
  pc <- paste0("prob_G", g)
  hist(assign_df[[pc]], breaks = 30, col = cols[g], border = "white",
       main = paste0("Group ", g, " Posterior Probability"),
       xlab = "Posterior Probability", xlim = c(0, 1))
  abline(v = 0.70, lty = 2, col = "gray40")
}
dev.off()

gc <- table(assign_df$trajectory_group)
png(file.path(DIRS$figures, "stage2_group_distribution.png"), width = 600, height = 500, res = 100)
barplot(gc, col = cols[names(gc)], border = "white",
        main = "Group Distribution", xlab = "Group", ylab = "N")
dev.off()

log_msg("   Figures saved")

log_msg("9. generating QA report...")
qa_lines <- c(
  "# Stage 2 QA: GBMTM trajectory modeling",
  "",
  paste0("**Time**: ", format(Sys.time())),
  paste0("**Selected model**: ng=", best_ng, " (d=2, scaling=0)"),
  "",
  "## Model comparison",
  "",
  "| ng | BIC | Entropy | APPA(mean) | MinGroup% | Converged |",
  "|----|-----|---------|------------|-----------|------|"
)
for (i in 1:nrow(compare_table)) {
  r <- compare_table[i, ]
  qa_lines <- c(qa_lines, sprintf("| %d | %.1f | %.4f | %.4f | %.1f%% | %s |",
    r$ng, r$BIC, r$entropy, r$APPA_mean, 100 * r$min_group_prop,
    ifelse(r$converged, "OK", "FAIL")))
}
qa_lines <- c(qa_lines, "",
  "## Per-group posterior diagnostics",
  sprintf("- Entropy: %.4f", compare_table$entropy[compare_table$ng == best_ng]),
  sprintf("- APPA(mean): %.4f", compare_table$APPA_mean[compare_table$ng == best_ng]),
  "",
  paste0("![Trajectory plot](../figures/stage2_trajectory_by_group.png)"),
  "",
  "> Stage 2 completed。")
writeLines(qa_lines, file.path(DIRS$qc_reports, "stage2_qc.md"))

log_msg("Stage 2 completed")
log_msg("============================================================")
