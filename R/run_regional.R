
initial_options <- commandArgs(trailingOnly = FALSE); file_arg <- initial_options[grep("--file=", initial_options)]; SCRIPTS_DIR <- if (length(file_arg) > 0) normalizePath(dirname(sub("--file=", "", file_arg[1]))) else normalizePath(".")
source(file.path(SCRIPTS_DIR, "config.R"))

bl <- read.csv(fp("baseline_5class.csv"), stringsAsFactors = FALSE)
bl$region <- factor(bl$region, levels = c("东北", "东部", "中部", "西部"))
tab <- table(bl$trajectory_group, bl$region, useNA = "no")
pct <- round(100 * prop.table(tab, 1), 1)
out <- as.data.frame.matrix(pct); out$n <- as.integer(rowSums(tab)); out$group <- paste0("G", rownames(out))
out <- out[, c("group", "n", "东北", "东部", "中部", "西部")]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
write.csv(out, file.path(OUT, "table_3.3_region.csv"), row.names = FALSE)
ct <- suppressWarnings(chisq.test(tab)); cat("chi-square p =", format.pval(ct$p.value, digits = 3), "\n")
png(file.path(OUT, "fig3_region_composition.png"), width = 2400, height = 1500, res = 300)
par(mar = c(4, 5, 2, 1))
m <- as.matrix(pct)
barplot(t(m), horiz = TRUE, col = c("#7FC8F8", "#5AA9E6", "#2E79B5", "#0B3C5D"),
        names.arg = paste0("G", rownames(pct), " (n=", as.integer(rowSums(tab)), ")"),
        xlab = "Share of class (%)", las = 1, cex.names = 0.9, border = NA)
legend("bottomright", c("Northeast", "East", "Central", "West"), fill = c("#7FC8F8", "#5AA9E6", "#2E79B5", "#0B3C5D"), bty = "n", cex = 0.8)
dev.off()
logmsg("written table_3.3_region.csv / fig3_region_composition.png")
