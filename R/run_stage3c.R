library(mice); library(dplyr)
initial_options <- commandArgs(trailingOnly = FALSE); file_arg <- initial_options[grep("--file=", initial_options)]; if (length(file_arg) > 0) { SCRIPTS_DIR <- normalizePath(dirname(sub("--file=", "", file_arg[1]))) } else { SCRIPTS_DIR <- normalizePath(dirname(sys.frame(1)$ofile)) }
source(file.path(SCRIPTS_DIR, "config.R"))

log_msg("============================================================")
log_msg("Stage 3C [Sandbox]: Missing Data Assessment + MICE Imputation")
log_msg("============================================================")

bl <- read.csv(file.path(DIRS$data_processed, "baseline_with_region.csv"), stringsAsFactors = FALSE)
anal_vars <- STAGE3$frozen_vars

miss_pct <- colMeans(is.na(bl[, anal_vars])) * 100
mice_vars <- names(miss_pct[miss_pct > 0])
log_msg(sprintf("Variables with missing: %d/%d", length(mice_vars), length(anal_vars)))
for (v in mice_vars) log_msg(sprintf("  %-30s %5.1f%%", label_var(v), miss_pct[v]))

df_imp <- bl[, c(anal_vars, "trajectory_group")]
df_imp$trajectory_group <- factor(df_imp$trajectory_group)

df_imp$region <- factor(df_imp$region, levels = c("东部","中部","西部","东北"))

meth <- rep("", ncol(df_imp))
names(meth) <- names(df_imp)
meth["trajectory_group"] <- ""
for (v in mice_vars) {
  if (v == "region") {
    meth[v] <- "polyreg"
  } else {
    uv <- unique(na.omit(df_imp[[v]]))
    if (all(uv %in% c(0, 1))) meth[v] <- "logreg"
    else meth[v] <- "pmm"
  }
}

log_msg(sprintf("Running MICE: m=%d, maxit=%d...", STAGE3$mice_m, STAGE3$mice_maxit))
imp <- mice(df_imp, m = STAGE3$mice_m, maxit = STAGE3$mice_maxit,
            method = meth, seed = STAGE3$mice_seed, print = FALSE)

log_msg("Post-processing and saving...")
for (m in 1:STAGE3$mice_m) {
  ds <- complete(imp, m)
  ds$trajectory_group <- as.integer(as.character(ds$trajectory_group))
  if ("region" %in% names(ds)) ds$region <- as.character(ds$region)
  for (v in mice_vars) {
    if (v != "region") {
      uv <- unique(na.omit(bl[[v]]))
      if (all(uv %in% c(0, 1))) ds[[v]] <- round(pmax(pmin(ds[[v]], 1), 0))
    }
  }
  write.csv(ds, file.path(DIRS$imputed, sprintf("imputed_m%02d.csv", m)), row.names = FALSE)
}
log_msg(sprintf("Saved %d imputed datasets (%d rows each)", STAGE3$mice_m, nrow(bl)))

log_msg("Generating density plot...")
pct_vars <- setdiff(mice_vars, "region")
png(file.path(DIRS$figures, "stage3c_imputation_density.png"), width = 1600, height = 1200, res = 130)
nr <- (length(pct_vars) + 3) %/% 4
if (nr > 0) {
  par(mfrow = c(nr, 4), mar = c(3, 3, 2, 1))
  for (v in pct_vars) {
    dens_obs <- density(bl[[v]], na.rm = TRUE)
    plot(dens_obs, main = sprintf("%s (%.1f%%)", label_var(v), miss_pct[v]), col = "steelblue", lwd = 2)
    for (mm in 1:min(3, STAGE3$mice_m)) {
      ds_m <- complete(imp, mm)
      lines(density(ds_m[[v]]), col = c("darkorange", "forestgreen", "firebrick")[mm], lwd = 1, lty = 2)
    }
  }
}
dev.off()

ds1 <- read.csv(file.path(DIRS$imputed, "imputed_m01.csv"))
n_miss <- sum(is.na(ds1))
log_msg(sprintf("Missing in imputed m01: %d (should be 0)", n_miss))

qa <- c("# Stage 3C QA: MICE Imputation [Sandbox]",
        "", paste0("**Time**: ", format(Sys.time())),
        paste0("**Method**: mice (pmm/logreg), m = ", STAGE3$mice_m, ", maxit = ", STAGE3$mice_maxit),
        paste0("**Variables imputed**: ", length(mice_vars)),
        paste0("**Remaining missing**: ", n_miss),
        "", "## Imputed Variables",
        paste(sprintf("- %s: %.1f%%", sapply(mice_vars, label_var), miss_pct[mice_vars]), collapse = "\n"),
        "", paste0("![Density plot](../figures/stage3c_imputation_density.png)"),
        "", "> Stage 3C [Sandbox] complete.")
writeLines(qa, file.path(DIRS$qc_reports, "stage3c_imputation_qc.md"))

log_msg("Stage 3C [Sandbox] complete")
