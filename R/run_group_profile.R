
initial_options <- commandArgs(trailingOnly = FALSE); file_arg <- initial_options[grep("--file=", initial_options)]; SCRIPTS_DIR <- if (length(file_arg) > 0) normalizePath(dirname(sub("--file=", "", file_arg[1]))) else normalizePath(".")
source(file.path(SCRIPTS_DIR, "config.R"))

m  <- readRDS(fp("gbmtm_ng5.rds"))
gr <- read.csv(fp("group_reordering.csv"), stringsAsFactors = FALSE)
asg<- read.csv(fp("final_trajectory_assignment.csv"), stringsAsFactors = FALSE)
n_by <- table(asg$trajectory_group)

post <- m$posterior; eps <- 1e-12; p <- pmax(pmin(post, 1 - eps), eps)
Ek <- 1 - sum(-rowSums(p * log(p))) / (nrow(p) * log(ncol(p)))

res <- data.frame()
for (g in sort(unique(gr$new_label))) {
  r  <- gr[gr$new_label == g, ]
  fit<- m$fitted[[r$group]]
  res <- rbind(res, data.frame(
    group = paste0("G", g), n = as.integer(n_by[as.character(g)]),
    pct = round(100 * as.integer(n_by[as.character(g)]) / sum(n_by), 1),
    PM_t0 = round(fit[1, 2], 1), PM_t2 = round(fit[2, 2], 1), PM_t4 = round(fit[3, 2], 1), PM_t7 = round(fit[4, 2], 1),
    PM_mean = round(mean(fit[, 2]), 1),
    Cog_t0 = round(fit[1, 1], 2), Cog_t2 = round(fit[2, 1], 2), Cog_t4 = round(fit[3, 1], 2), Cog_t7 = round(fit[4, 1], 2),
    Cog_change = round(fit[4, 1] - fit[1, 1], 2),
    Cog_range = round(max(fit[, 1]) - min(fit[, 1]), 2),
    risk = round(r$risk_score, 3), stringsAsFactors = FALSE))
}
res$APPA <- round(as.numeric(m$appa)[match(gr$group[order(gr$new_label)], seq_along(m$appa))], 3)
res$Entropy <- round(Ek, 3)
res$mean_pp <- round(tapply(asg$posterior_prob, asg$trajectory_group, mean)[as.character(1:5)], 3)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
write.csv(res, file.path(OUT, "table_3.3_group_profile.csv"), row.names = FALSE)

png(file.path(OUT, "fig2_trajectory_5class.png"), width = 3000, height = 1400, res = 300)
par(mfrow = c(1, 2), mar = c(4.2, 4.4, 2, 1))
tp <- 1:4; tl <- c("2011", "2013", "2015", "2018")
cols <- c("#1B4965", "#3E92CC", "#5FA8D3", "#C1666B", "#7D4E57")
pchs <- c(15, 16, 17, 18, 19)
for (k in 1:2) {
  y <- if (k == 1) sapply(gr$group[order(gr$new_label)], function(i) m$fitted[[i]][, 1]) else
                    sapply(gr$group[order(gr$new_label)], function(i) m$fitted[[i]][, 2])
  plot(tp, y[, 1], type = "n", xaxt = "n", ylim = range(y),
       xlab = "Wave", ylab = if (k == 1) "Cognition (raw score)" else "PM2.5 (ug/m3)", box = TRUE)
  axis(1, at = tp, labels = tl)
  for (g in 1:5) lines(tp, y[, g], col = cols[g], lwd = 2.2, type = "o", pch = pchs[g])
  legend(if (k == 1) "bottomleft" else "topright", paste0("G", 1:5), col = cols, pch = pchs, lwd = 2, cex = 0.85, bty = "n")
}
dev.off()
logmsg("written table_3.3_group_profile.csv / fig2_trajectory_5class.png")
