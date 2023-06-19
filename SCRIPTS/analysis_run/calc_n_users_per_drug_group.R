# DESCRIPTION: Calculate the number of users per drug group
# AUTHOR: Julia Romanowska
# DATE CREATED: 2023-01-30
# DATE MODIFIED: 2023-05-22

# SETUP ----
library(tidyverse)
library(here)
library(data.table)
library(survival)

time_lag <- 0
prescriptions_exposure <- 8

dataset_atc2level_file <- paste0("dataset_ready_for_analysis_exposure",
																 prescriptions_exposure,".rds")
users_per_drug_file <- paste0(
	"n_users_per_ATC_code2_exposure", prescriptions_exposure,
	"_time-lag", time_lag, "yrs.txt"
)
users_per_drug_no_strata_file <- paste0(
	"n_users_per_ATC_code2_no_sex_exposure", prescriptions_exposure,
	"_time-lag", time_lag, "yrs.txt"
)

# READ DATA ----
data_atc2level <- readRDS(
	here("DATA", dataset_atc2level_file)
)
data_atc2level <- data_atc2level[age >= 24]

dim(data_atc2level)
names(data_atc2level)

## how many cases, how many controls?
data_atc2level[,.N, by = park_yn]

## PREPARE FOR TIME LAG ----
if(time_lag > 0){
	data_atc2level[
		, time_risk := time_risk - time_lag*365.25
	]
}

# CHECK USAGE ----
all_atc_codes <- str_subset(
	string = colnames(data_atc2level),
	pattern = "^[[:upper:][:digit:][:digit:]]"
)
all_atc_codes

other_colnames <- setdiff(
	colnames(data_atc2level),
	all_atc_codes
)
other_colnames

all_atc_codes <- str_subset(all_atc_codes, pattern = "^V0", negate = TRUE)

n_users_per_atc_group <- map(all_atc_codes, function(cur_atc_code){
	cur_columns <- c(other_colnames, cur_atc_code)
	cur_data <- as_tibble(
		data_atc2level[, ..cur_columns]
		) %>%
		rename(mtid = all_of(cur_atc_code)) %>%
		# only those who had the first prescription _after_ 2004
		filter(time_risk > 0) %>%
		# fix problems with negative time
		mutate(mtid = if_else(mtid < 0, 0.5, mtid))
	
	print(cur_data %>% count(park_yn))
	print(
		cur_data %>%
		filter(!is.na(mtid)) %>%
		count(park_yn)
	)
	
	all_n_users_non_pd <- cur_data %>%
		filter(!is.na(mtid) & !park_yn) %>%
		count(park_yn)
	
	# cur_data %>%
	# 	filter(!is.na(mtid) & park_yn & mtid >= time_risk) %>% glimpse()
	## THAT'S WHAT I WANT TO REMOVE: PD patients that began using the drug
	##    too late
	all_n_users_pd <- cur_data %>%
		filter(!is.na(mtid) & park_yn & mtid < time_risk) %>%
		count(park_yn)
	
	sex_strat_n_users_non_pd <- cur_data %>%
		filter(!is.na(mtid) & !park_yn) %>%
		count(sex, park_yn)
	
	sex_strat_n_users_pd <- cur_data %>%
		filter(!is.na(mtid) & park_yn & mtid < time_risk) %>%
		count(sex, park_yn)
	
	# combine:
	return(
		list(
			all = 
				bind_rows(all_n_users_pd, all_n_users_non_pd) %>%
				rename(n_users = n) %>%
				add_column(ATC_code = cur_atc_code),
			sex_strat = 
				bind_rows(sex_strat_n_users_non_pd, sex_strat_n_users_pd) %>%
				rename(n_users = n) %>%
				add_column(ATC_code = cur_atc_code)
		)
	)
})

n_users_no_strat <- map(
	n_users_per_atc_group,
	~ pluck(.x, "all")
) %>%
	bind_rows()

write_delim(
	n_users_no_strat,
	file = here("DATA", users_per_drug_no_strata_file),
	delim = "\t"
)
n_users_no_strat

n_users_sex_strat <- map(
	n_users_per_atc_group,
	~ pluck(.x, "sex_strat")
) %>%
	bind_rows()
n_users_sex_strat

write_delim(
	n_users_sex_strat,
	file = here("DATA", users_per_drug_file),
	delim = "\t"
)
