library(haven); library(dplyr)
initial_options <- commandArgs(trailingOnly = FALSE); file_arg <- initial_options[grep("--file=", initial_options)]; if (length(file_arg) > 0) { SCRIPTS_DIR <- normalizePath(dirname(sub("--file=", "", file_arg[1]))) } else { SCRIPTS_DIR <- normalizePath(dirname(sys.frame(1)$ofile)) }
source(file.path(SCRIPTS_DIR, "config.R"))

log_msg("============================================================")
log_msg("Stage 3A: Baseline Profiling")
log_msg("============================================================")

log_msg("Loading CHARLS 2011...")
df_raw <- read_dta(PATHS$charls_main)
health_2011 <- read_dta(PATHS$health_2011)
for (col in names(df_raw))
  if (inherits(df_raw[[col]], "haven_labelled"))
    df_raw[[col]] <- haven::zap_labels(df_raw[[col]])
for (col in names(health_2011))
  if (inherits(health_2011[[col]], "haven_labelled"))
    health_2011[[col]] <- haven::zap_labels(health_2011[[col]])
df_raw$ID    <- sprintf("%012.0f", as.numeric(df_raw$ID))
df_raw$ID_w1 <- as.character(as.numeric(df_raw$ID_w1))
health_2011$ID <- as.character(as.numeric(health_2011$ID))

log_msg("Mapping variables...")
df_bl <- data.frame(ID = df_raw$ID, ID_w1 = df_raw$ID_w1)

df_bl$sex <- ifelse(df_raw$ragender == 1, 1, 0)
edu_raw  <- df_raw$raeduc_c
df_bl$education <- dplyr::case_when(
  edu_raw <= 3 ~ 1, edu_raw == 4 ~ 2, edu_raw == 5 ~ 3,
  edu_raw %in% 6:7 ~ 4, edu_raw >= 8 ~ 5, TRUE ~ NA_real_)
df_bl$marital_status <- ifelse(df_raw$r1mstat %in% c(1, 3), 1, 0)
df_bl$residence <- ifelse(df_raw$h1rural == 1, 0, 1)
df_bl$smoking   <- df_raw$r1smokev
df_bl$drinking  <- df_raw$r1drinkev
df_bl$bmi       <- df_raw$r1mbmi
df_bl$self_rated_health <- 6 - df_raw$r1shlt
df_bl$adl_score <- df_raw$r1adla_c
df_bl$depression_cesd <- df_raw$r1cesd10

df_bl$hypertension <- df_raw$r1hibpe;  df_bl$diabetes <- df_raw$r1diabe
df_bl$dyslipidemia <- df_raw$r1dyslipe; df_bl$heart_disease <- df_raw$r1hearte
df_bl$stroke <- df_raw$r1stroke;       df_bl$lung_disease <- df_raw$r1lunge
df_bl$asthma <- df_raw$r1asthmae;      df_bl$arthritis <- df_raw$r1arthre
df_bl$cancer <- df_raw$r1cancre;       df_bl$liver_disease <- df_raw$r1livere
df_bl$kidney_disease <- df_raw$r1kidneye
df_bl$digestive_disease <- df_raw$r1digeste
df_bl$mobility_impairment <- df_raw$r1mobilsev
df_bl$employment <- ifelse(df_raw$r1work == 1, 1, 0)
df_bl$age <- df_raw$r1agey

log_msg("Health module...")
df_hl <- health_2011[, c("ID", "da049", "da050")]
names(df_hl) <- c("ID_w1_hl", "sleep_night", "nap")
df_hl$sleep_night <- as.numeric(df_hl$sleep_night)
df_hl$nap <- as.numeric(df_hl$nap)
df_bl <- left_join(df_bl, df_hl, by = c("ID_w1" = "ID_w1_hl"))

log_msg("Housing characteristics (solid fuel)...")
df_housing <- read_dta(file.path(dirname(PATHS$health_2011), "housing_characteristics.dta"))
for (col in names(df_housing))
  if (inherits(df_housing[[col]], "haven_labelled"))
    df_housing[[col]] <- haven::zap_labels(df_housing[[col]])
df_housing$solid_fuel <- ifelse(df_housing$i022 %in% c(1, 6), 1,
                                ifelse(df_housing$i022 %in% c(2, 3, 4, 5), 0, NA_real_))
df_housing$householdID <- as.character(as.numeric(df_housing$householdID))
df_bl$householdID_w1 <- as.character(as.numeric(df_raw$householdID_w1))
df_bl <- left_join(df_bl, df_housing[, c("householdID", "solid_fuel")],
                    by = c("householdID_w1" = "householdID"))
df_bl$householdID_w1 <- NULL
log_msg(sprintf("Solid fuel: %d users, %d clean, %d NA",
                sum(df_bl$solid_fuel == 1, na.rm = TRUE),
                sum(df_bl$solid_fuel == 0, na.rm = TRUE),
                sum(is.na(df_bl$solid_fuel))))

df_bl$bmi[df_bl$bmi < 10 | df_bl$bmi > 60] <- NA
log_msg(sprintf("BMI cleaned: %d outliers -> NA", sum(df_bl$bmi < 10 | df_bl$bmi > 60, na.rm = TRUE)))

df_bl$metabolic_burden       <- rowSums(cbind(df_bl$hypertension, df_bl$diabetes, df_bl$dyslipidemia), na.rm = TRUE)
df_bl$cardio_cerebrovascular <- rowSums(cbind(df_bl$heart_disease, df_bl$stroke), na.rm = TRUE)
df_bl$respiratory_burden     <- rowSums(cbind(df_bl$lung_disease, df_bl$asthma), na.rm = TRUE)

log_msg("Merging trajectory assignment...")
assign_df <- read.csv(file.path(DIRS$models, "final_trajectory_assignment.csv"))
assign_df$ID <- sprintf("%012.0f", as.numeric(as.character(assign_df$ID)))
df_bl <- left_join(df_bl, assign_df[, c("ID", "trajectory_group", "posterior_prob")], by = "ID")
df_bl <- df_bl[!is.na(df_bl$trajectory_group), ]
anal_vars <- STAGE3$frozen_vars
anal_vars <- setdiff(anal_vars, "region")
log_msg(sprintf("Merged: %d individuals x %d columns", nrow(df_bl), ncol(df_bl)))

log_msg("Generating Table 1...")
groups <- sort(unique(df_bl$trajectory_group))
t1 <- data.frame(variable = sapply(anal_vars, label_var), stringsAsFactors = FALSE)
for (g in groups) {
  col_name <- paste0("G", g)
  t1[[col_name]] <- sapply(anal_vars, function(v) {
    vals <- df_bl[[v]][df_bl$trajectory_group == g]
    sprintf("%.2f (%.2f)", mean(vals, na.rm = TRUE), sd(vals, na.rm = TRUE))
  })
}
p_raw <- sapply(anal_vars, function(v) {
  vals <- df_bl[[v]]; grp <- factor(df_bl$trajectory_group)
  if (v %in% c("age", "bmi", "sleep_night", "nap")) {
    oneway.test(vals ~ grp, var.equal = FALSE)$p.value
  } else if (v %in% c("education", "self_rated_health", "adl_score",
                       "depression_cesd", "mobility_impairment",
                       "metabolic_burden", "cardio_cerebrovascular",
                       "respiratory_burden")) {
    kruskal.test(vals ~ grp)$p.value
  } else {
    chisq.test(table(vals, grp))$p.value
  }
})
p_adj <- p.adjust(p_raw, method = "bonferroni")
t1$P_value <- sapply(p_adj, function(p) {
  if (p < 0.001) "<0.001" else sprintf("%.3f", p)
})
write.csv(t1, file.path(DIRS$tables, "table1_baseline_characteristics.csv"), row.names = FALSE)
log_msg(sprintf("Table 1: %d vars, P<0.05 (Bonferroni-adjusted): %d",
                nrow(t1), sum(p_adj < 0.05, na.rm = TRUE)))

write.csv(df_bl, file.path(DIRS$data_processed, "baseline_with_trajectory.csv"), row.names = FALSE)
log_msg("Baseline saved")
log_msg("Stage 3A complete")
