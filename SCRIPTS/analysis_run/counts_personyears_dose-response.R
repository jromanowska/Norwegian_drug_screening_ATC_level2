# DESCRIPTION: Create counts of person-years for dose-response tables
# AUTHOR: Julia Romanowska
# DATE CREATED: 2023-06-06
# DATE LAST MODIFIED: 

# SETUP ----
library(data.table)
library(tidyverse)
library(here)
library(gt)
library(gtsummary)
library(survival)

# data with status date and covariates
dataset_atc2level_file <- paste0(
	"dataset_ready_for_analysis_exposure2.rds"
)
# data with "dose" dates
dataset_atc2level_dose_dates_file <- "time_dependent_drug_count_usage_ATC_level2.rds"

# READ DATA ----
data_atc2level <- readRDS(
	here("DATA", dataset_atc2level_file)
)

dosage_dates <- readRDS(
	here("DATA", dataset_atc2level_dose_dates_file)
)

# FILTER -----------
# remove all persons who were 23 or younger in 2004
data_atc2level <- data_atc2level[age >= 24]

data_atc2level[
	, age0 := (age * 365.25) + 365.25/2
]

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

personal_data <- data_atc2level[
	, ..other_colnames
]

# COUNT ----
threshold_labels <- paste0("Q", 1:5)
start_date <- lubridate::ymd("2005-01-01")

count_personyears_atc2 <- map(
	atc2names,
	function(atc){
		cat("\n === CURRENT ATC GROUP: ", atc, " ===\n")
		cur_dosage_dates <- dosage_dates[[atc]]
		
	# first, I need to recreate the data used for cox:
	cur_dosage_data <- dosage_dates[[atc]]
	if(is.null(cur_dosage_data)){
		return(NULL)
	}
	# transform dosage data so that it includes 4 columns with dates for
	#   each dosage crossing
	cur_dosage_data[, prescr_date := lubridate::ymd(
		paste(prescr_year, prescr_mnth, "15", sep = "-")
	)][
		, prescr_time := as.numeric(prescr_date - ..start_date)
	]
	cur_dosage_data_wide <- dcast(
		cur_dosage_data,
		indiv_ID ~ cur_threshold,
		value.var = "prescr_time"
	)
	
	# I need to add Q1 explicitly if it's not there - this happens when first
	#   step has the same number of prescriptions as the definition of exposure
	if(!("Q1" %in% names(cur_dosage_data_wide))){
		cur_dosage_data_wide[
			, Q1 := Q2
		]
	}
	
	cur_data_input <- merge(
		personal_data,
		cur_dosage_data_wide,
		by = "indiv_ID",
		all.x = TRUE
	)
	cur_data_input <- cur_data_input[
		time_risk > 0,
	]
	cur_data_input[
		, time_risk := time_risk + age0
	][
		, c(threshold_labels) := lapply(
			.SD,
			function(x){
				if_else(x < 0, 0.5, x)
			}
		),
		.SDcols = threshold_labels
	][
		, Q2 := if_else(
			!is.na(Q2) & Q1 == Q2,
			Q2 + 0.05,
			Q2
		)
	][
		, Q3 := if_else(
			!is.na(Q3) & Q2 == Q3,
			Q3 + 0.05,
			Q3
		)
	][
		, Q4 := if_else(
			!is.na(Q4) & Q3 == Q4,
			Q4 + 0.05,
			Q4
		)
	][
		, Q5 := if_else(
			!is.na(Q5) & Q3 == Q5,
			Q5 + 0.05,
			Q5
		)
	][
		, c(threshold_labels) := lapply(
			.SD,
			function(x) x + age0
		),
		.SDcols = threshold_labels
	]
	
	cur_data <- tmerge(
		data1 = cur_data_input,
		data2 = cur_data_input,
		id = indiv_ID,
		park = event(time_risk, park_yn),
		trt = tdc(Q1),
		trt = tdc(Q2),
		trt = tdc(Q3),
		trt = tdc(Q4),
		trt = tdc(Q5)
	)
	cur_data <- tmerge(
		cur_data,
		cur_data,
		id = indiv_ID,
		trt2 = cumtdc(tstart)
	)
	cur_data <- cur_data %>% 
		mutate(tstart = if_else(is.na(Q1), 
														tstart + age0, tstart)) %>% 
		mutate(tstart = if_else(tstart == 0, 
														tstart + age0, tstart))
	# now the data can be used to count:
	
	## from Julia:
	person_years <- as_tibble(cur_data) %>%
		mutate(py = (tstop - tstart) / 365.25) %>% 
		select(indiv_ID, py, trt2) %>% 
		pivot_wider(
			names_from = "trt2",
			values_from = "py",
			names_prefix = "trt_q") %>% 
		setnafill(., fill = 0) %>%
		tbl_summary(include = -c("indiv_ID"),
								statistic = all_continuous() ~ "{sum}",
								digits = all_continuous() ~ 0)
	#person_years
	
	cases <- as_tibble(cur_data) %>%
		select(indiv_ID, park, trt2) %>%
		mutate(park = if_else(
			park, 1, 0
		)) %>%
		pivot_wider(
			names_from = "trt2",
			values_from = "park",
			names_prefix = "trt_q") %>% 
		setnafill(., fill = 0) %>% 
		tbl_summary(include = -c("indiv_ID"),
								statistic = all_dichotomous() ~ "{n}")
	#cases
	
	cur_personyears_case_perdose <- tbl_merge(
		tbls = list(person_years, cases),
		tab_spanner = c("**Personyears**", "**Cases**")
		)
	cur_personyears_case_perdose
	
	return(as_tibble(cur_personyears_case_perdose))
}
)
names(count_personyears_atc2) <- atc2names

# save results ----
saveRDS(
	count_personyears_atc2,
	file = here("DATA", "person_years_and_cases_per_dose_ATC2level.rds")
)
