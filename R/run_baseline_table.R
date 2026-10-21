
initial_options <- commandArgs(trailingOnly = FALSE); file_arg <- initial_options[grep("--file=", initial_options)]; SCRIPTS_DIR <- if (length(file_arg) > 0) normalizePath(dirname(sub("--file=", "", file_arg[1]))) else normalizePath(".")
source(file.path(SCRIPTS_DIR, "config.R"))

bl <- read.csv(fp("baseline_5class.csv"), stringsAsFactors = FALSE)
WELCH <- c("age","bmi","sleep_night","nap")
KRUSK <- c("education","self_rated_health","adl_score","depression_cesd",
           "mobility_impairment","metabolic_burden","cardio_cerebrovascular","respiratory_burden")
BIN <- setdiff(VARS, c(WELCH, KRUSK))
grp <- factor(bl$trajectory_group)
p <- sapply(VARS, function(v) {
  x <- bl[[v]]
  if (v %in% WELCH) oneway.test(x ~ grp, var.equal = FALSE)$p.value
  else if (v %in% KRUSK) kruskal.test(x ~ grp)$p.value
  else chisq.test(table(x, grp))$p.value
})
padj <- p.adjust(p, method = "bonferroni")
out <- data.frame(variable = LAB[VARS], stringsAsFactors = FALSE)
for (g in 1:5) {
  s <- bl[bl$trajectory_group == g, ]
  out[[paste0("G", g)]] <- sapply(VARS, function(v) {
    x <- s[[v]]
    if (v %in% BIN) sprintf("%d (%.1f)", sum(x == 1, na.rm = TRUE), 100 * mean(x == 1, na.rm = TRUE))
    else sprintf("%.2f (%.2f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
  })
}
out$P_adj <- ifelse(padj < 0.001, "<0.001", sprintf("%.3f", padj))
out$P_raw <- signif(p, 3)
ns <- sapply(1:5, function(g) sum(bl$trajectory_group == g))
attr(out, "N") <- ns
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
write.csv(out, file.path(OUT, "table_3.4_baseline_5class.csv"), row.names = FALSE)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
write.csv(data.frame(group = paste0("G", 1:5), N = ns), file.path(OUT, "table_3.4_N.csv"), row.names = FALSE)
cat("significant (Bonferroni<0.05):", sum(padj < 0.05), "/", length(VARS), "\n")
logmsg("written table_3.4_baseline_5class.csv")
