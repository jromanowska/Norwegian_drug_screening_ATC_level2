# DESCRIPTION: Prepare data for dose-response analyses;
#   the number of prescriptions will be used as a proxy for dose;
#   this is because we're doing analysis on level2, which groups different
#   drugs into one, so the actual dose value from the prescriptions would not
#   give meaningful results.
# AUTHOR: Julia Romanowska
# DATE CREATED: 2023-06-01
# DATE LAST MODIFIED: 2023-06-02

# SETUP ----
library(data.table)
library(tidyverse)
library(DrugRegScreening)
library(here)
library(dbplyr)

exposure_threshold <- 2

# READ DATA ----
dataset_atc2level_file <- paste0(
	"dataset_ready_for_analysis_exposure", exposure_threshold, ".rds"
)
data_atc2level <- readRDS(
	here("DATA", dataset_atc2level_file)
)
# remove all persons who were 23 or younger in 2004
data_atc2level <- data_atc2level[age >= 24]
data_atc2level
dim(data_atc2level)
data_atc2level[,.N, by = park_yn]

# get all the names of drug groups
atc2names <- str_subset(
	string = colnames(data_atc2level),
	pattern = "^[[:upper:][:digit:][:digit:]]"
)
atc2names

other_colnames <- setdiff(
	colnames(data_atc2level),
	atc2names
)
other_colnames

# database with all data
my_db <- getDBcon()

setDT(all_prescriptions <- tbl(my_db, "Prescription_data") %>%
		select(indiv_ID:no_doses, ATC_code, starts_with("refusjon"),
					 starts_with("DDD")) %>%
		collect())
gc()

all_prescriptions[
	, ATC_code2 := str_sub(ATC_code, start = 1, end = 3)
]
gc()
setkey(all_prescriptions, indiv_ID, ATC_code2)
gc()

all_prescriptions

# PREPARE TIME-DEPT.DATA ----
quartile_vector <- c(Q1 = 0.25, Q2 = 0.5, Q3 = 0.75, Q4 = 0.9)

time_dept_dose_usage <- map(
	atc2names,
	function(cur_atc_group){
		cat("\n ====== CURRENT ATC GROUP: ", cur_atc_group, "===== \n")
		
		cur_prescriptions_dt <- all_prescriptions[
			ATC_code2 == cur_atc_group
		]
		setkey(cur_prescriptions_dt, indiv_ID, prescr_year)
		
		## FILTER PRESCRIPTIONS ----
		# leave only those that were taken before censoring
		status_date_per_indiv <- data_atc2level[
			, .(indiv_ID, status_date, park_yn)
		]
		cur_prescriptions_dt[
			, prescr_date := lubridate::ymd(
				paste(prescr_year, prescr_mnth, "15", sep = "-")
			)
		]
		filtered_prescriptions_dt <- cur_prescriptions_dt[
			status_date_per_indiv,
			on = .(indiv_ID = indiv_ID, prescr_date < status_date),
			nomatch = NULL
		]
		#filtered_prescriptions_dt
		
		prescriptions_per_indiv_dt <- filtered_prescriptions_dt[
			, .N,
			by = .(indiv_ID, prescr_year, prescr_mnth)
		]
		#prescriptions_per_indiv_dt
		prescriptions_per_indiv_dt[
			, dummy_N := 1
		][
			, cumsum_N := cumsum(dummy_N),
			by = indiv_ID
		]
		max_n_prescr_per_indiv_dt <- prescriptions_per_indiv_dt[
			, max_N := max(cumsum_N),
			by = indiv_ID
		][
			cumsum_N == max_N,
		]
		
		all_exposed <- max_n_prescr_per_indiv_dt[
			cumsum_N >= exposure_threshold,
			.(indiv_ID, prescr_year, prescr_mnth, cumsum_N)
		]
		cur_thresholds <- quantile(all_exposed$cumsum_N, quartile_vector)
		cat("   current thresholds: ", as.numeric(cur_thresholds), "\n")
		# check whether several breaks are the same - too few prescriptions
		if(length(unique(as.numeric(cur_thresholds))) <
			 length(cur_thresholds)){
			return(NULL)
		}

		thresholded_prescriptions_per_indiv_dt <- prescriptions_per_indiv_dt[
			cumsum_N >= exposure_threshold,
		]
		thresholded_prescriptions_per_indiv_dt[
			, cur_threshold := santoku::chop(cumsum_N, cur_thresholds)
		]
		thresholded_prescriptions_per_indiv_dt <- thresholded_prescriptions_per_indiv_dt[
			, .SD[1], by = .(indiv_ID, cur_threshold)
		]
		#thresholded_prescriptions_per_indiv_dt
		distrib_before <- summary(thresholded_prescriptions_per_indiv_dt$cur_threshold)
		thresholded_prescriptions_per_indiv_dt[
			, cur_threshold := fct_relabel(
				cur_threshold,
				function(x){
					case_when(
						x == "25%" ~ "Q2",
						x == "50%" ~ "Q3",
						x == "75%" ~ "Q4",
						x == "90%" ~ "Q5",
						is.na(x) ~ NA_character_,
						TRUE ~ "Q1"
					)
				}
			)
		]
		distrib_after <- summary(thresholded_prescriptions_per_indiv_dt$cur_threshold)
		stopifnot(identical(as.numeric(distrib_before), as.numeric(distrib_after)))
		return(
			thresholded_prescriptions_per_indiv_dt[
				, .(indiv_ID, cur_threshold, prescr_year, prescr_mnth, cumsum_N, max_N)
			]
		)
	}
)

names(time_dept_dose_usage) <- atc2names

saveRDS(
	time_dept_dose_usage,
	here("DATA", "time_dependent_drug_count_usage_ATC_level2.rds")
)

closeDBcon(my_db)
