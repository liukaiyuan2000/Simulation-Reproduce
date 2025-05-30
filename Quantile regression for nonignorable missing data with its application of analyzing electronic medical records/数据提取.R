# 加载必要的包
library(dplyr)
library(readr)
library(lubridate)
# 设置数据路径
data_path <- "D:/Desktop/mimic-iii-clinical-database-1.4/"

# 读取需要的表格
patients <- read_csv(paste0(data_path, "PATIENTS.csv"))
admissions <- read_csv(paste0(data_path, "ADMISSIONS.csv"))
icustays <- read_csv(paste0(data_path, "ICUSTAYS.csv"))
patients <- read_csv(paste0(data_path, "PATIENTS.csv"))
patients$AGE <- year(ymd(patients$DOD)) - year(ymd(patients$DOB))

# 合并表格并筛选符合条件的患者
filtered_data <- patients %>%
  # 合并住院信息
  inner_join(admissions, by = "SUBJECT_ID") %>%
  # 合并ICU住院信息
  inner_join(icustays, by = c("SUBJECT_ID", "HADM_ID")) %>%
  
  # 筛选种族为黑人或白人
  filter(ETHNICITY %in% c("BLACK/AFRICAN AMERICAN", "WHITE")) %>%
  
  # 筛选有保险的患者（政府或私人保险）
  filter(INSURANCE %in% c("Government", "Private")) %>%
  
  # 筛选有完整记录的患者
  filter(
    !is.na(GENDER),
    !is.na(DOB),  # 出生日期，用于计算年龄
    !is.na(MARITAL_STATUS),
    !is.na(LANGUAGE)
  ) %>%
  
  # 筛选婚姻状态为已婚或单身
  filter(MARITAL_STATUS %in% c("MARRIED", "SINGLE")) %>%
  
  # 处理语言信息：英语或非英语
  mutate(
    FIRST_LANGUAGE = case_when(
      grepl("ENGLISH", LANGUAGE, ignore.case = TRUE) ~ "English",
      TRUE ~ "non-English"
    )
  ) %>%
  
  # 确保年龄合理（过滤掉异常值）
  filter(AGE >= 1 & AGE <= 150) %>%  # 根据研究需求调整
  
  # 选择需要的列
  select(
    SUBJECT_ID,
    GENDER,
    AGE,
    MARITAL_STATUS,
    FIRST_LANGUAGE,
    INSURANCE,
    ETHNICITY
  ) %>%
  
  # 去除重复患者（每人只保留一条记录）
  distinct(SUBJECT_ID, .keep_all = TRUE)

# 检查最终数据量
message(paste("筛选后患者数量:", nrow(filtered_data)))



# chartevents <- read_csv(paste0(data_path, "CHARTEVENTS.csv"), 
#                         col_types = cols_only(
#                           SUBJECT_ID = col_integer(),
#                           HADM_ID = col_integer(),
#                           ICUSTAY_ID = col_integer(),
#                           ITEMID = col_integer(),
#                           CHARTTIME = col_character(),
#                           VALUENUM = col_double(),
#                           VALUEUOM = col_character()
#                         ))












