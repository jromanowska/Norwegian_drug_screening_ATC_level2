# DESCRIPTION: check some descriptive stats of the sample
# AUTHOR: Julia Romanowska
# DATE CREATED: 2022-10-05
# DATE MODIFIED: 2023-02-07

# SETUP ----
library(tidyverse)
library(here)
library(data.table)
library(dbplyr)
library(RSQLite)
library(DBI)
library(gtsummary)

dataset_atc2level_file <- "dataset_ready_for_analysis.rds"


# READ DATA ----
data_atc2level <- readRDS(
	here("DATA", dataset_atc2level_file)
)
data_atc2level <- data_atc2level[age >= 24]

dim(data_atc2level)
names(data_atc2level)

# CHECKING STATS ----
## descriptive table ----

## how many cases, how many controls? ----
data_atc2level[,.N, by = park_yn]

## mean age at PD diagnosis ----
# diagnosis date = 'first_N04_date'
# 'age' column is the age on 1.1.2005

# extract the necessary info
data_atc2level_personal <- data_atc2level[
	, .(indiv_ID, park_yn, first_N04_date,
			time_to_onset, foedselsaar, sex, age, edu_level,
			censor_event, cens_date, status_date, time_risk)
]

age_at_diagn <- as_tibble(data_atc2level_personal[
	park_yn == TRUE,
])
age_at_diagn %>%
	select(indiv_ID, park_yn, first_N04_date, age)

age_at_diagn <- age_at_diagn %>%
	mutate(age_at_onset = age + (year(first_N04_date) - 2005))
age_at_diagn %>%
	select(indiv_ID, park_yn, first_N04_date, age, age_at_onset)

# mean age:
age_at_diagn %>%
	summarise(mean_age_at_onset = mean(age_at_onset))

## how many male patients? ----
sex_freq_all <- data_atc2level_personal[, .N, by = .(sex)]
sex_freq_all[
	, prcnt := N/sum(N) * 100
]
sex_freq_all

sex_freq <- data_atc2level_personal[, .N, by = .(park_yn, sex)]
sex_freq[
	, prcnt := N/sum(N) * 100
]
sex_freq

#only PD patients:
sex_freq_PD <- data_atc2level_personal[park_yn == TRUE, .N, by = sex]
sex_freq_PD[
	, prcnt := N/sum(N) * 100
]
sex_freq_PD

## nice descriptive table? ----
descr_data <- data_atc2level_personal[
	, .(park_yn, age, sex, edu_level)
]
descr_data[
	, sex_female := if_else(
		sex == 2, TRUE, FALSE
	)
]
descr_data[
	, edu_level := if_else(
		edu_level == 0, NA_real_, edu_level
	)
][
	, edu_level := factor(
		edu_level,
		levels = c(NA, 1:4),
		labels = c("primary school", "highschool", "college/university (short)", "college/university (long)")
	)
]

descr_table <- tbl_summary(
	descr_data[, -"sex"],
	by = "park_yn",
	label = list(
		age ~ "Age at start of follow-up",
		sex_female ~ "Female sex",
		edu_level ~ "Education level"
	)
)
descr_table %>%
	as_gt() %>%
	gt::gtsave(filename = here("RESULTS", "descriptive_table.tex"))

## how many had received levodopa, how many MAO-B? ----
# I will fetch this info from the DB, based on the indiv_ID from age_at_diagn
db_con <- dbConnect(RSQLite::SQLite(), "../../DrugScreeningDB.sqlite")
# Get all prescriptions for the patients

if(FALSE){
# to make it faster: create index first!
	# need to run it only once!
RSQLite::dbExecute(
	db_con,
	"CREATE INDEX indiv_idx ON Prescription_data( indiv_ID )"
)
RSQLite::dbExecute(
	db_con,
	"CREATE INDEX ATC_idx ON Prescription_data( ATC_code )"
)
}

# CHECKING:
dbGetQuery(
	db_con,
	"SELECT indiv_ID, prescr_year, prescr_mnth, ATC_code
	FROM Prescription_data
	WHERE ATC_code LIKE 'N04BA%'
	LIMIT 10"
)

all_prescr_mao_b_levo <- as_tibble(
	dbGetQuery(
		db_con,
		"SELECT indiv_ID, prescr_year, prescr_mnth, ATC_code
				FROM Prescription_data
				WHERE ATC_code LIKE 'N04BD%' OR ATC_code LIKE 'N04BA%'"
	)
)
all_prescr_mao_b_levo %>% count(ATC_code)

count_prescr_patients <- age_at_diagn %>%
	left_join(
		## --- THIS DIDN'T WORK! ---
		# tbl(db_con, "Prescription_data") %>%
		# 	select(indiv_ID:prescr_mnth, ATC_code, refusjonICD:refusjonICPC) %>%
		# 	filter(ATC_code %like% "^NO4BD" | ATC_code %like% "^N04BA"),
		## ---
		all_prescr_mao_b_levo,
		by = c("indiv_ID" = "indiv_ID")
	) %>%
	mutate(levo_or_mao_b = case_when(
		ATC_code %like% "^N04BD" ~ "mao_b",
		ATC_code %like% "^N04BA" ~ "levo",
		TRUE ~"other"
	)) %>%
	count(indiv_ID, levo_or_mao_b)

count_prescr_patients %>%
	select(-n) %>%
	count(levo_or_mao_b) %>%
	mutate(prcnt = n/nrow(age_at_diagn)*100)

dbDisconnect(db_con)
