library(haven)
library(readxl)
library(dplyr)
library(tidyr)

initial_options <- commandArgs(trailingOnly = FALSE); file_arg <- initial_options[grep("--file=", initial_options)]; if (length(file_arg) > 0) { SCRIPTS_DIR <- normalizePath(dirname(sub("--file=", "", file_arg[1]))) } else { SCRIPTS_DIR <- normalizePath(dirname(sys.frame(1)$ofile)) }
source(file.path(SCRIPTS_DIR, "config.R"))

log_msg("============================================================")
log_msg("Stage 1: Data preprocessing start")
log_msg("============================================================")

log_msg("1. reading CHARLS main data...")
df_raw <- read_dta(PATHS$charls_main)
log_msg(sprintf("   raw: %d rows x %d cols, %d individuals",
                nrow(df_raw), ncol(df_raw), n_distinct(df_raw$ID)))

log_msg("2. computing cognitive score...")

clean_draw <- function(x) {
  x <- as.character(x)
  x[x == "1. Yes" | x == "1.Yes"] <- "1"
  x[x == "0. No"  | x == "0.No"]  <- "0"
  x[grepl("^\\.", x)] <- NA_character_
  as.numeric(x)
}

for (wave in names(STAGE1$cog_components)) {
  comps <- STAGE1$cog_components[[wave]]
  draw_col <- grep("draw", comps, value = TRUE)
  if (length(draw_col) > 0 && draw_col %in% names(df_raw)) {
    df_raw[[draw_col]] <- clean_draw(df_raw[[draw_col]])
  }
  for (comp in comps) {
    if (comp %in% names(df_raw)) {
      df_raw[[comp]] <- as.numeric(df_raw[[comp]])
    }
  }
  total_col <- paste0(wave, "_total")
  df_raw[[total_col]] <- rowSums(df_raw[, comps, drop = FALSE], na.rm = FALSE)
  log_msg(sprintf("   %s: n=%d, mean=%.2f, sd=%.2f",
    wave, sum(!is.na(df_raw[[total_col]])),
    mean(df_raw[[total_col]], na.rm = TRUE),
    sd(df_raw[[total_col]], na.rm = TRUE)))
}

log_msg("3. converting to long format...")
total_cols <- paste0(names(STAGE1$cog_components), "_total")
df_wide <- df_raw[, c("ID", "ID_w1", total_cols)]

df_long <- df_wide %>%
  pivot_longer(cols = all_of(total_cols),
               names_to = "year_code", values_to = "cog_score") %>%
  mutate(
    year_num = as.integer(gsub("_total", "", gsub("r", "", year_code))),
    year = as.integer(STAGE1$year_mapping[as.character(year_num)])
  ) %>%
  select(ID, ID_w1, year, cog_score) %>%
  arrange(ID, year)
log_msg(sprintf("   long table: %d rows", nrow(df_long)))

log_msg("4. extracting age data...")
age_cols <- c("r1agey", "r2agey", "r3agey", "r4agey")
age_map <- c(r1agey = 2011, r2agey = 2013, r3agey = 2015, r4agey = 2018)
df_age <- df_raw[, c("ID", age_cols)] %>%
  pivot_longer(cols = all_of(age_cols),
               names_to = "year_code", values_to = "age") %>%
  mutate(year = as.integer(age_map[year_code])) %>%
  select(ID, year, age) %>%
  arrange(ID, year)
log_msg(sprintf("   age rows: %d", nrow(df_age)))

log_msg("5. address matching...")

match_address <- function(year, health_path, psu_path = NULL) {
  df_indiv <- read_dta(health_path)
  df_core <- distinct(df_indiv[, c("ID", "communityID")])

  if (!is.null(psu_path) && file.exists(psu_path)) {
    df_psu <- read_dta(psu_path)
    psu_map <- distinct(df_psu[, c("communityID", "province", "city",
                                    "province_eng", "city_eng")])
  } else {
    df_psu_2013 <- read_dta(PATHS$psu_2013)
    psu_map <- distinct(df_psu_2013[, c("communityID", "province", "city",
                                         "province_eng", "city_eng")])
  }
  df_result <- left_join(df_core, psu_map, by = "communityID") %>%
    select(-communityID) %>%
    mutate(year = year)
  matched <- sum(!is.na(df_result$province))
  log_msg(sprintf("   %d: total=%d, matched=%d (%.1f%%)",
    year, nrow(df_result), matched, 100*matched/nrow(df_result)))
  df_result
}

df_addr <- bind_rows(
  match_address(2011, PATHS$health_2011, PATHS$psu_2011),
  match_address(2013, PATHS$health_2013, PATHS$psu_2013),
  match_address(2015, PATHS$health_2015, NULL),
  match_address(2018, PATHS$health_2018, NULL)
)

log_msg("6. merging cognition, address, age...")
df_cog_2011 <- filter(df_long, year == 2011)
df_addr_2011 <- filter(df_addr, year == 2011) %>% rename(ID_w1 = ID)
df_merged_2011 <- left_join(df_cog_2011, df_addr_2011, by = c("ID_w1", "year"))

df_cog_other <- filter(df_long, year != 2011)
df_addr_other <- filter(df_addr, year != 2011)
df_merged_other <- left_join(df_cog_other, df_addr_other, by = c("ID", "year"))

df_merged <- bind_rows(df_merged_2011, df_merged_other) %>%
  arrange(ID, year) %>%
  left_join(df_age, by = c("ID", "year"))

log_msg(sprintf("   merged: %d rows, %d individuals", nrow(df_merged), n_distinct(df_merged$ID)))

log_msg("7. city/province name mapping...")
df_merged <- df_merged %>%
  mutate(
    city_cn     = CITY_MAP[city_eng],
    province_cn = PROVINCE_MAP[province_eng],
    city        = ifelse(is.na(city_cn), as.character(city), city_cn),
    province    = ifelse(is.na(province_cn), as.character(province), province_cn)
  ) %>%
  select(-city_cn, -province_cn)
log_msg(sprintf("   mapping completed"))

log_msg("8. matching PM2.5 pollution data...")
df_poll <- read_excel(PATHS$pollution)
names(df_poll) <- c("province_city", "year", "O3", "PM25", "PM10", "CO", "NO2", "SO2")
df_poll[df_poll == "No Data"] <- NA_character_
for (col in c("PM25", "PM10", "O3", "CO", "NO2", "SO2")) {
  df_poll[[col]] <- as.numeric(df_poll[[col]])
}

parse_region <- function(name) {
  if (name == "中华人民共和国") {
    return(c("national", "全国"))
  } else if (grepl("省$", name)) {
    if (nchar(name) == 3) {
      return(c("province", name))
    } else {
      return(c("city", name))
    }
  } else if (grepl("自治区$", name)) {
    return(c("province", name))
  } else if (grepl("市$", name)) {
    return(c("city", name))
  } else if (grepl("州$|盟$|地区$", name)) {
    return(c("city", name))
  } else {
    return(c("unknown", name))
  }
}

parsed <- t(sapply(df_poll$province_city, parse_region))
df_poll$region_level <- parsed[, 1]
df_poll$region_name  <- parsed[, 2]

unknown_regions <- unique(df_poll$region_name[df_poll$region_level == "unknown"])
if (length(unknown_regions) > 0) {
  log_msg(sprintf("   warning: %d regions unclassifiable: %s",
    length(unknown_regions), paste(unknown_regions, collapse = ", ")))
}

match_pm25 <- function(row_year, row_city, row_province, df_pol) {
  city_match <- df_pol[df_pol$year == row_year &
                       df_pol$region_level == "city" &
                       df_pol$region_name == row_city, ]
  if (nrow(city_match) > 0) {
    return(c(city_match$PM25[1], "城市匹配"))
  }
  prov_match <- df_pol[df_pol$year == row_year &
                       df_pol$region_level %in% c("province", "city") &
                       df_pol$region_name == row_province, ]
  if (nrow(prov_match) > 0) {
    return(c(prov_match$PM25[1], "省份匹配"))
  }
  return(c(NA_character_, "未匹配"))
}

n_rows <- nrow(df_merged)
pm_results <- matrix(NA_character_, nrow = n_rows, ncol = 2)
for (i in seq_len(n_rows)) {
  if (i %% 10000 == 0) log_msg(sprintf("   PM2.5 match progress: %d/%d", i, n_rows))
  res <- match_pm25(df_merged$year[i], df_merged$city[i],
                    df_merged$province[i], df_poll)
  pm_results[i, ] <- res
}
df_merged$PM25 <- as.numeric(pm_results[, 1])
df_merged$match_source <- pm_results[, 2]
log_msg(sprintf("   PM2.5 matching completed: %s", paste(table(pm_results[, 2]), collapse = ", ")))

log_msg("9. applying inclusion/exclusion criteria...")
initial_n <- n_distinct(df_merged$ID)
log_msg(sprintf("   before screening: %d individuals", initial_n))

df_age_ok <- df_merged %>%
  group_by(ID) %>%
  filter(all(age >= STAGE1$inclusion$min_age)) %>%
  ungroup()
n_age <- n_distinct(df_age_ok$ID)
log_msg(sprintf("   age >=%d: %d (excluded %d)",
  STAGE1$inclusion$min_age, n_age, initial_n - n_age))

df_selected <- df_age_ok %>%
  group_by(ID) %>%
  filter(sum(!is.na(cog_score)) >= STAGE1$inclusion$min_waves_cog &
         sum(!is.na(PM25)) >= STAGE1$inclusion$min_waves_pm25) %>%
  ungroup()
n_cog_pm <- n_distinct(df_selected$ID)
log_msg(sprintf("    cognition waves>=%d & PM2.5 waves>=%d: %d (excluded %d)",
  STAGE1$inclusion$min_waves_cog, STAGE1$inclusion$min_waves_pm25,
  n_cog_pm, n_age - n_cog_pm))

df_selected <- df_selected %>%
  group_by(ID) %>%
  filter(all(match_source == "城市匹配")) %>%
  ungroup()
final_n <- n_distinct(df_selected$ID)
log_msg(sprintf("    city-match only: %d (excluded %d)",
  final_n, n_cog_pm - final_n))
log_msg(sprintf("    final: %d individuals, %d rows", final_n, nrow(df_selected)))

log_msg("10. time encoding...")
df_selected <- df_selected %>%
  mutate(time = as.numeric(STAGE1$time_encoding[as.character(year)]))
stopifnot(!any(is.na(df_selected$time)))

log_msg("11. saving output...")
df_out <- df_selected %>%
  select(ID, year, time, cog_score, PM25) %>%
  arrange(ID, year)

output_path <- file.path(DIRS$data_processed, STAGE1$output_file)
write.csv(df_out, output_path, row.names = FALSE)
log_msg(sprintf("   output: %s", output_path))
log_msg(sprintf("   dimensions: %d rows x 5 cols", nrow(df_out)))
log_msg(sprintf("   individuals: %d", n_distinct(df_out$ID)))
log_msg(sprintf("   years: %s", paste(names(table(df_out$year)), collapse = ", ")))

log_msg("12. city consistency check...")
city_stay <- df_selected %>%
  group_by(ID) %>%
  summarise(
    n_waves = n(),
    n_cities = n_distinct(city),
    all_same = n_distinct(city) == 1,
    .groups = "drop"
  )
n_same_city <- sum(city_stay$all_same)
n_total_ids <- nrow(city_stay)
log_msg(sprintf("   consistent city across all waves: %d/%d (%.1f%%)",
                n_same_city, n_total_ids, 100 * n_same_city / n_total_ids))
city_counts <- table(city_stay$n_cities)
c1 <- ifelse(is.na(city_counts["1"]), 0, city_counts["1"])
c2 <- ifelse(is.na(city_counts["2"]), 0, city_counts["2"])
c3 <- ifelse(is.na(city_counts["3"]), 0, city_counts["3"])
c4 <- ifelse(is.na(city_counts["4"]), 0, city_counts["4"])
log_msg(sprintf("   cities per individual: 1=%d (%.1f%%), 2=%d, 3=%d, 4=%d",
                c1, 100*c1/n_total_ids, c2, c3, c4))

log_msg("Stage 1 completed")
log_msg("============================================================")
