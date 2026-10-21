set.seed(2026)

SCRIPT_DIR <- tryCatch({
  initial_options <- commandArgs(trailingOnly = FALSE)
  file_arg <- initial_options[grep("--file=", initial_options)]
  if (length(file_arg) > 0) {
    normalizePath(dirname(sub("--file=", "", file_arg[1])))
  } else {
    normalizePath(dirname(sys.frame(1)$ofile))
  }
}, error = function(e) {
  normalizePath(file.path(getwd(), "sandbox", "scripts"))
})
SANDBOX_ROOT <- dirname(SCRIPT_DIR)

.find_project_root <- function(start) {
  d <- normalizePath(start, winslash = "/", mustWork = FALSE)
  for (i in seq_len(8)) {
    if (dir.exists(file.path(d, "data"))) return(d)
    p <- dirname(d)
    if (identical(p, d)) break
    d <- p
  }
  normalizePath(file.path(start, "..", "..", "..", ".."), winslash = "/", mustWork = FALSE)
}
PROJECT_ROOT <- .find_project_root(SCRIPT_DIR)

DATA_DIR <- file.path(PROJECT_ROOT, "data", "charls")
FALLBACK_CHARLS_DIR <- Sys.getenv("CHARLS_DIR", unset = "data/charls")
FALLBACK_POLLUTION  <- Sys.getenv("POLLUTION_FILE", unset = "data/pollution.xlsx")

.dp <- function(rel) {
  p <- file.path(DATA_DIR, rel)
  if (file.exists(p)) p else file.path(FALLBACK_CHARLS_DIR, rel)
}
.dp_pollution <- function() {
  p <- file.path(DATA_DIR, "空气污染数据", "Combined_Pollutants_with_Units.xlsx")
  if (file.exists(p)) p else FALLBACK_POLLUTION
}

PATHS <- list(
  charls_main = .dp("2011-2018/H_CHARLS_D_Data.dta"),
  health_2011 = .dp("2011/2011health_status_and_functioning.dta"),
  health_2013 = .dp("2013/2013Health_Status_and_Functioning.dta"),
  health_2015 = .dp("2015/2015Health_Status_and_Functioning.dta"),
  health_2018 = .dp("2018/Health_Status_and_Functioning.dta"),
  psu_2011    = .dp("2011/2011PSU.dta"),
  psu_2013    = .dp("2013/2013PSU.dta"),
  pollution   = .dp_pollution()
)

DIRS <- list(
  data_processed = file.path(SANDBOX_ROOT, "output", "data"),
  tables         = file.path(SANDBOX_ROOT, "output", "tables"),
  figures        = file.path(SANDBOX_ROOT, "output", "figures"),
  models         = file.path(SANDBOX_ROOT, "output", "models"),
  logs           = file.path(SANDBOX_ROOT, "output", "logs"),
  qc_reports     = file.path(SANDBOX_ROOT, "output", "qc_reports"),
  imputed        = file.path(SANDBOX_ROOT, "output", "imputed")
)
for (d in DIRS) dir.create(d, showWarnings = FALSE, recursive = TRUE)

STAGE1 <- list(
  cog_components = list(
    r1 = c("r1orient", "r1tr20", "r1draw", "r1ser7"),
    r2 = c("r2orient", "r2tr20", "r2draw", "r2ser7"),
    r3 = c("r3orient", "r3tr20", "r3draw", "r3ser7"),
    r4 = c("r4orient", "r4tr20", "r4draw", "r4ser7")
  ),
  year_mapping = c("1" = 2011, "2" = 2013, "3" = 2015, "4" = 2018),
  time_encoding = c("2011" = 0, "2013" = 2, "2015" = 4, "2018" = 7),
  inclusion = list(min_age = 45, min_waves_cog = 3, min_waves_pm25 = 3),
  output_file = "trajectory_long.csv"
)

STAGE2 <- list(
  d = 2,
  scaling = 0,
  pruning = TRUE,
  ng_range = c(2, 3, 4, 5, 6),
  selection = list(min_entropy = 0.70, min_appa = 0.70, min_group_prop = 0.10),
  n_cores = max(1, parallel::detectCores() - 1)
)

STAGE3 <- list(
  reference_group = 1,
  frozen_vars = c(
    "age", "sex", "education", "marital_status", "residence",
    "smoking", "drinking", "sleep_night", "nap", "employment", "solid_fuel",
    "bmi", "self_rated_health", "adl_score", "depression_cesd",
    "arthritis", "cancer", "liver_disease", "kidney_disease", "digestive_disease",
    "mobility_impairment", "metabolic_burden", "cardio_cerebrovascular", "respiratory_burden"
  ),
  mice_m = 10,
  mice_maxit = 5,
  mice_seed = 2026
)

VAR_LABELS <- c(
  age = "Age", sex = "Sex", education = "Education",
  marital_status = "Marital status", residence = "Residence",
  employment = "Employment",
  smoking = "Smoking", drinking = "Drinking",
  sleep_night = "Sleep night", nap = "Nap",
  bmi = "BMI", self_rated_health = "Self-rated health",
  adl_score = "ADL score", depression_cesd = "Depression (CESD-10)",
  mobility_impairment = "Mobility impairment",
  arthritis = "Arthritis", cancer = "Cancer",
  liver_disease = "Liver disease", kidney_disease = "Kidney disease",
  digestive_disease = "Digestive disease",
  metabolic_burden = "Metabolic burden",
  cardio_cerebrovascular = "Cardio-cerebrovascular",
  respiratory_burden = "Respiratory burden",
  hypertension = "Hypertension", diabetes = "Diabetes",
  dyslipidemia = "Dyslipidemia",
  heart_disease = "Heart disease", stroke = "Stroke",
  lung_disease = "Lung disease", asthma = "Asthma",
  trajectory_group = "Trajectory group",
  posterior_prob = "Posterior probability",
  cog_score = "Cognitive score", PM25 = "PM2.5",
  solid_fuel = "Solid fuel use"
)

label_var <- function(v) {
  ifelse(v %in% names(VAR_LABELS), VAR_LABELS[v], v)
}

PROVINCE_MAP <- c(
  Yunnan = "云南省", Fujian = "福建省", Qinghai = "青海省",
  Sichuan = "四川省", Hebei = "河北省", Jiangxi = "江西省",
  Xinjiang = "新疆维吾尔自治区", Beijing = "北京市",
  Neimenggu = "内蒙古自治区", Jiangsu = "江苏省",
  Chongqing = "重庆市", Gansu = "甘肃省",
  Heilongjiang = "黑龙江省", Guangdong = "广东省",
  Liaoning = "辽宁省", Shanxi = "山西省", Shanghai = "上海市",
  Tianjin = "天津市", Zhejiang = "浙江省", Jilin = "吉林省",
  Guangxi = "广西壮族自治区", Anhui = "安徽省", Hubei = "湖北省",
  Shaanxi = "陕西省", Shandong = "山东省", Henan = "河南省",
  Hunan = "湖南省", Guizhou = "贵州省"
)

CITY_MAP <- c(
  Baoshan = "保山市", Chuxiong = "楚雄彝族自治州", Zhaotong = "昭通市",
  Kunming = "昆明市", Lincang = "临沧市", Lijiang = "丽江市",
  Zhangzhou = "漳州市", Putian = "莆田市", Fuzhou = "福州市",
  Ningde = "宁德市", Haidong = "海东市", Guangan = "广安市",
  Chengdu = "成都市", Liangshan = "凉山彝族自治州",
  Nanchong = "南充市", Yibin = "宜宾市", Ziyang = "资阳市",
  Neijiang = "内江市", Ganzi = "甘孜藏族自治州",
  Mianyang = "绵阳市", Meishan = "眉山市",
  Baoding = "保定市", Shijiazhuang = "石家庄市",
  Cangzhou = "沧州市", Chengde = "承德市",
  Jingdezhen = "景德镇市", Nanchang = "南昌市",
  Yichun = "宜春市", Shangrao = "上饶市", Jiujiang = "九江市",
  Jian = "吉安市", Ganzhou = "赣州市",
  Akesu = "阿克苏地区", Beijing = "北京市",
  Xilingol = "锡林郭勒盟", Hohhot = "呼和浩特市",
  Hinggan = "兴安盟", Chifeng = "赤峰市", Hulunbuir = "呼伦贝尔市",
  Suzhou = "苏州市", Xuzhou = "徐州市", Taizhou = "泰州市",
  Yancheng = "盐城市", Suqian = "宿迁市", Yangzhou = "扬州市",
  Lianyungang = "连云港市", Chongqing = "重庆市",
  Lanzhou = "兰州市", Dingxi = "定西市", Pingliang = "平凉市",
  Zhangye = "张掖市", Jixi = "鸡西市", Qiqihar = "齐齐哈尔市",
  Harbin = "哈尔滨市", Jiamusi = "佳木斯市",
  Foshan = "佛山市", Shenzhen = "深圳市", Chaozhou = "潮州市",
  Guangzhou = "广州市", Maoming = "茂名市",
  Qingyuan = "清远市", Jiangmen = "江门市",
  Benxi = "本溪市", Anshan = "鞍山市", Dalian = "大连市",
  Chaoyang = "朝阳市", Jinzhou = "锦州市",
  Yangquan = "阳泉市", Xinzhou = "忻州市", Yuncheng = "运城市",
  Linfen = "临汾市", Shanghai = "上海市", Tianjin = "天津市",
  Huzhou = "湖州市", Ningbo = "宁波市", Hangzhou = "杭州市",
  Lishui = "丽水市", Jiaxing = "嘉兴市",
  Siping = "四平市", Jilin = "吉林市",
  Guilin = "桂林市", Hechi = "河池市", Nanning = "南宁市",
  Yulin = "玉林市", Bozhou = "亳州市", Fuyang = "阜阳市",
  Huainan = "淮南市", Anqing = "安庆市", Liuan = "六安市",
  Chaohu = "巢湖市", Xiangfan = "襄阳市",
  Enshi = "恩施土家族苗族自治州", Huanggang = "黄冈市",
  Jingmen = "荆门市", Weinan = "渭南市", Baoji = "宝鸡市",
  Hanzhong = "汉中市", Binzhou = "滨州市", Qingdao = "青岛市",
  Jinan = "济南市", Zaozhuang = "枣庄市", Linyi = "临沂市",
  Liaocheng = "聊城市", Weihai = "威海市", Weifang = "潍坊市",
  Dezhou = "德州市", Anyang = "安阳市", Luoyang = "洛阳市",
  Zhoukou = "周口市", Zhengzhou = "郑州市", Puyang = "濮阳市",
  Pingdingshan = "平顶山市", Jiaozuo = "焦作市",
  Xinyang = "信阳市", Shaoyang = "邵阳市", Yueyang = "岳阳市",
  Changsha = "长沙市", Yiyang = "益阳市", Loudi = "娄底市",
  Changde = "常德市", Qiannan = "黔南布依族苗族自治州",
  Qiandongnan = "黔东南苗族侗族自治州"
)

LOG_PATH <- file.path(DIRS$logs, "pipeline.log")
log_msg <- function(msg) {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  line <- paste(timestamp, msg)
  cat(line, "\n")
  write(line, LOG_PATH, append = TRUE)
  flush.console()
}

STAGE3$reference_group <- 2

STAGE3$frozen_vars <- c("region", STAGE3$frozen_vars)
VAR_LABELS["region"] <- "Region"

PROVINCE_TO_REGION <- list(
  east = c(
    "北京市", "天津市", "河北省", "上海市", "江苏省",
    "浙江省", "福建省", "山东省", "广东省", "海南省"
  ),
  central = c(
    "山西省", "安徽省", "江西省", "河南省", "湖北省", "湖南省"
  ),
  west = c(
    "内蒙古自治区", "广西壮族自治区", "重庆市", "四川省",
    "贵州省", "云南省", "西藏自治区", "陕西省", "甘肃省",
    "青海省", "宁夏回族自治区", "新疆维吾尔自治区"
  ),
  northeast = c(
    "辽宁省", "吉林省", "黑龙江省"
  )
)

get_region <- function(province_chinese) {
  if (province_chinese %in% PROVINCE_TO_REGION$east) return("东部")
  if (province_chinese %in% PROVINCE_TO_REGION$central) return("中部")
  if (province_chinese %in% PROVINCE_TO_REGION$west) return("西部")
  if (province_chinese %in% PROVINCE_TO_REGION$northeast) return("东北")
  return(NA_character_)
}

log_msg("Config loaded: REF=2, region in MICE only.")
log_msg(paste("Project root:", PROJECT_ROOT))
log_msg(paste("output:", SANDBOX_ROOT))
log_msg("Seed: 2026")

LAB <- c(sex="Sex (male)", residence="Urban residence", smoking="Current smoking",
         drinking="Current drinking", marital_status="Married/partnered",
         employment="Employed", solid_fuel="Solid fuel use", arthritis="Arthritis",
         cancer="Cancer", liver_disease="Liver disease", kidney_disease="Kidney disease",
         digestive_disease="Digestive disease", age="Age (years)", bmi="BMI (kg/m2)",
         self_rated_health="Self-rated health (1-5)", adl_score="ADL score (0-4)",
         depression_cesd="Depression (CESD-10)", mobility_impairment="Mobility impairment (0-7)",
         metabolic_burden="Metabolic burden (0-3)", cardio_cerebrovascular="Cardio-cerebrovascular (0-2)",
         respiratory_burden="Respiratory burden (0-2)", education="Education (level 1-5)",
         sleep_night="Nocturnal sleep (h)", nap="Nap (min)")
VARS <- names(LAB)

IN <- DIRS$data_processed
OUT <- DIRS$tables
logmsg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
fp <- function(f) {
  for (dd in c(DIRS$data_processed, DIRS$models, DIRS$tables)) {
    p <- file.path(dd, f)
    if (file.exists(p)) return(p)
  }
  file.path(IN, f)
}
