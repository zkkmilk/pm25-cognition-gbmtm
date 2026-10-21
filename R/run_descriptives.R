
initial_options <- commandArgs(trailingOnly = FALSE); file_arg <- initial_options[grep("--file=", initial_options)]; SCRIPTS_DIR <- if (length(file_arg) > 0) normalizePath(dirname(sub("--file=", "", file_arg[1]))) else normalizePath(".")
source(file.path(SCRIPTS_DIR, "config.R"))

long <- read.csv(fp("trajectory_long.csv"), stringsAsFactors = FALSE)
bl   <- read.csv(fp("baseline_5class.csv"), stringsAsFactors = FALSE)
logmsg("3.1 descriptives: n=", nrow(bl))

f <- function(x, dig = 1) sprintf(paste0("%.", dig, "f (%.", dig, "f)"), mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
rows <- list(
  data.frame(item = "N", value = as.character(nrow(bl))),
  data.frame(item = "Age, mean (SD)", value = f(bl$age)),
  data.frame(item = "Male, n (%)", value = sprintf("%d (%.1f)", sum(bl$sex == 1), 100 * mean(bl$sex == 1, na.rm = TRUE))),
  data.frame(item = "Rural, n (%)", value = sprintf("%d (%.1f)", sum(bl$residence == 0, na.rm = TRUE), 100 * mean(bl$residence == 0, na.rm = TRUE)))
)
edu <- table(factor(bl$education, levels = 1:5))
edu_lab <- c("No formal / did not finish primary", "Primary", "Junior high", "High school / vocational", "College or above")
for (i in 1:5) rows[[length(rows) + 1]] <- data.frame(item = paste0("Education: ", edu_lab[i]),
                                                       value = sprintf("%d (%.1f)", edu[i], 100 * edu[i] / sum(edu)))
tab <- do.call(rbind, rows)

trend <- aggregate(cbind(PM25, cog_score) ~ time, data = long, FUN = function(x) c(m = mean(x, na.rm = TRUE), s = sd(x, na.rm = TRUE)))
pm <- aggregate(PM25 ~ time, long, function(x) c(m = mean(x, na.rm = TRUE), s = sd(x, na.rm = TRUE)))
cg <- aggregate(cog_score ~ time, long, function(x) c(m = mean(x, na.rm = TRUE), s = sd(x, na.rm = TRUE)))
tr <- data.frame(time = pm$time, PM25_mean = round(pm$PM25[, "m"], 1), PM25_sd = round(pm$PM25[, "s"], 1),
                 Cog_mean = round(cg$cog_score[, "m"], 2), Cog_sd = round(cg$cog_score[, "s"], 2))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
write.csv(tab, file.path(OUT, "table_3.1_sample.csv"), row.names = FALSE)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
write.csv(tr,  file.path(OUT, "table_3.1_trends.csv"), row.names = FALSE)
logmsg("written table_3.1_sample.csv / table_3.1_trends.csv")
