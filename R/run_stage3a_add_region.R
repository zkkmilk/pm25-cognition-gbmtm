library(haven)
library(dplyr)

initial_options <- commandArgs(trailingOnly = FALSE); file_arg <- initial_options[grep("--file=", initial_options)]; if (length(file_arg) > 0) { SCRIPTS_DIR <- normalizePath(dirname(sub("--file=", "", file_arg[1]))) } else { SCRIPTS_DIR <- normalizePath(dirname(sys.frame(1)$ofile)) }
source(file.path(SCRIPTS_DIR, "config.R"))

log_msg("step region")

baseline_path <- file.path(DIRS$data_processed, "baseline_with_trajectory.csv")
baseline <- read.csv(baseline_path, stringsAsFactors = FALSE)
log_msg(sprintf("读取 baseline: %d 人, %d 列", nrow(baseline), ncol(baseline)))

health_2011 <- read_dta(PATHS$health_2011)
health_ids <- health_2011 %>%
  select(ID, communityID) %>%
  distinct() %>%
  mutate(ID = as.character(trimws(ID)),
         communityID = as.character(communityID))

log_msg(sprintf("健康模块: %d 唯一ID, %d 唯一communityID",
                n_distinct(health_ids$ID), n_distinct(health_ids$communityID)))
log_msg(sprintf("健康模块 ID 样本: %s", paste(head(health_ids$ID, 5), collapse = ", ")))

baseline <- baseline %>% mutate(
  ID_w1_padded = ifelse(nchar(ID_w1) == 10, paste0("0", ID_w1), ID_w1)
)
log_msg(sprintf("baseline ID_w1(原) 样本: %s", paste(head(baseline$ID_w1, 5), collapse = ", ")))
log_msg(sprintf("baseline ID_w1(补零) 样本: %s", paste(head(baseline$ID_w1_padded, 5), collapse = ", ")))

id_to_comm <- baseline %>%
  select(ID, ID_w1_padded) %>%
  left_join(health_ids, by = c("ID_w1_padded" = "ID")) %>%
  select(ID, communityID)

n_comm_match <- sum(!is.na(id_to_comm$communityID))
log_msg(sprintf("ID_w1→communityID 匹配: %d / %d (%.1f%%)",
                n_comm_match, nrow(id_to_comm), 100 * n_comm_match / nrow(id_to_comm)))

psu11 <- read_dta(PATHS$psu_2011)

psu_prov <- psu11 %>%
  select(communityID, province_eng) %>%
  mutate(
    communityID = as.character(communityID),
    province_chinese = PROVINCE_MAP[province_eng],
    region = sapply(province_chinese, get_region)
  )

n_psu <- n_distinct(psu_prov$communityID)
n_id_comm <- sum(!is.na(id_to_comm$communityID))
n_match_comm <- length(intersect(
  unique(id_to_comm$communityID[!is.na(id_to_comm$communityID)]),
  unique(psu_prov$communityID)))
log_msg(sprintf("PSU communityID: %d, ID→comm 有效: %d, 交集: %d",
                n_psu, n_id_comm, n_match_comm))

missing_prov <- psu_prov %>% filter(is.na(province_chinese) | is.na(region))
if (nrow(missing_prov) > 0) {
  log_msg(sprintf("警告: %d PSU 记录省份/region 未匹配", nrow(missing_prov)))
}

comm_to_region <- id_to_comm %>%
  left_join(psu_prov %>% select(communityID, province_chinese, region), by = "communityID") %>%
  select(ID, province_chinese, region)

unmatched <- id_to_comm %>%
  filter(!communityID %in% psu_prov$communityID) %>%
  distinct(communityID)
n_unmatched_comm <- nrow(unmatched)
n_unmatched_ind <- sum(is.na(comm_to_region$region))
if (n_unmatched_comm > 0) {
  log_msg(sprintf("警告: %d 个 communityID 未匹配到 PSU, 涉及 %d 人",
    n_unmatched_comm, n_unmatched_ind))
  log_msg(sprintf("  未匹配 communityID: %s",
    paste(head(unmatched$communityID, 10), collapse = ", ")))
}

baseline_new <- baseline %>%
  left_join(comm_to_region %>% select(ID, province_chinese, region), by = "ID")

n_merged <- sum(!is.na(baseline_new$region))
log_msg(sprintf("region 匹配成功: %d / %d (%.1f%%)",
                n_merged, nrow(baseline_new), 100 * n_merged / nrow(baseline_new)))

log_msg("Region 分布:")
log_msg(sprintf("  东部(参照): %d (%.1f%%)",
  sum(baseline_new$region == "东部", na.rm = TRUE),
  100 * sum(baseline_new$region == "东部", na.rm = TRUE) / n_merged))
log_msg(sprintf("  中部: %d (%.1f%%)",
  sum(baseline_new$region == "中部", na.rm = TRUE),
  100 * sum(baseline_new$region == "中部", na.rm = TRUE) / n_merged))
log_msg(sprintf("  西部: %d (%.1f%%)",
  sum(baseline_new$region == "西部", na.rm = TRUE),
  100 * sum(baseline_new$region == "西部", na.rm = TRUE) / n_merged))
log_msg(sprintf("  东北: %d (%.1f%%)",
  sum(baseline_new$region == "东北", na.rm = TRUE),
  100 * sum(baseline_new$region == "东北", na.rm = TRUE) / n_merged))

baseline_new$ID_w1_padded <- NULL
baseline_new$province_chinese <- NULL

out_path <- file.path(DIRS$data_processed, "baseline_with_region.csv")
write.csv(baseline_new, out_path, row.names = FALSE)
log_msg(sprintf("已保存: %s (%d 人, %d 列)", out_path, nrow(baseline_new), ncol(baseline_new)))
log_msg("step region")
