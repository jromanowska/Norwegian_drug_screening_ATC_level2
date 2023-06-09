# DESCRIPTION: Dose-response analyses; looking at number of prescriptions per
#   ATC level2 code and per individual
# AUTHOR: Julia Romanowska
# DATE CREATED: 2023-06-02
# DATE LAST MODIFIED: 2023-06-05

# SETUP ----
library(data.table)
library(tidyverse)
library(here)
library(survival)
library(future)
library(future.apply)

exposure_threshold <- 2
n_parallel_processes <- 10

dirty_res_file <- here(
	"RESULTS", "dirty_results_atc_group2_time-lag0yrs_dose-response.rds"
)
tidy_res_file <- here(
	"RESULTS", "nice_results_atc_group2_time-lag0yrs_dose-response"
)

# data with status date and covariates
dataset_atc2level_file <- paste0(
	"dataset_ready_for_analysis_exposure", exposure_threshold, ".rds"
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
dim(data_atc2level)
data_atc2level[,.N, by = park_yn]

# categorise age into 5-year intervals
data_atc2level[
	, age_group := cut(
			age,
			breaks = c(seq(25, 85, by = 5), 115),
			right = FALSE
		)
]
data_atc2level

# calculate age (in days!) at the time of start of follow up
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

# ANALYSE DATA ----
if(!file.exists(dirty_res_file)){
	# prepare for parallelization
	options(future.globals.maxSize = 3 *10^9 * 1024^2)
	plan(multisession, workers = n_parallel_processes)
	
	start.time <- Sys.time()
	start.time
	
	personal_data <- data_atc2level[
		, ..other_colnames
	]
	
	cox_res_atc2 <- future_lapply(
		atc2names,
		function(atc){
			threshold_labels <- paste0("Q", 1:5)
			start_date <- lubridate::ymd("2005-01-01")
			
			cat("\n === CURRENT ATC GROUP: ", atc, " ===\n")
			
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
			#cur_data_input
			
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
			# head(cur_data, 20)
			cur_data <- tmerge(
				cur_data,
				cur_data,
				id = indiv_ID,
				trt2 = cumtdc(tstart)
			)
			
			# checking with cumevent
			if(FALSE){
			cur_data_input_long <- merge(
				personal_data,
				cur_dosage_data,
				by = "indiv_ID",
				all.x = TRUE
			)
			cur_data_input_long <- cur_data_input_long[
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
			cur_data_input
			
			cur_data_cumevent <- tmerge(
				cur_data
			)
			}
			
			cur_data <- cur_data %>% 
				mutate(tstart = if_else(is.na(Q1), 
																tstart + age0, tstart)) %>% 
				mutate(tstart = if_else(tstart == 0, 
																tstart + age0, tstart))
			
			cur_cox_res_all <- coxph(
				Surv(tstart, tstop, park) ~ 
					as.factor(trt2) + as.factor(edu_level) + as.factor(sex),
				data = cur_data
			)
			# cur_cox_res_all
			cur_survfit <- survfit(
				Surv(tstart, tstop, park) ~ 
					as.factor(trt2),
				data = cur_data
			)
			#survminer::ggsurvplot(cur_survfit, title = atc)
			
			cur_cox_res_women <- coxph(
				Surv(tstart, tstop, park) ~ 
					as.factor(trt2) + as.factor(edu_level),
				data = cur_data %>% filter(sex == 2)
			)
			cur_survfit_women <- survfit(
				Surv(tstart, tstop, park) ~ 
					as.factor(trt2),
				data = cur_data %>% filter(sex == 2)
			)
			
			cur_cox_res_men <- coxph(
				Surv(tstart, tstop, park) ~ 
					as.factor(trt2) + as.factor(edu_level),
				data = cur_data %>% filter(sex == 1)
			)
			cur_survfit_men <- survfit(
				Surv(tstart, tstop, park) ~ 
					as.factor(trt2),
				data = cur_data %>% filter(sex == 1)
			)
			
			return(
				list(
					cox_all = cur_cox_res_all,
					cox_women = cur_cox_res_women,
					cox_men = cur_cox_res_men,
					surv_all = cur_survfit,
					surv_women = cur_survfit_women,
					surv_men = cur_survfit_men
				)
			)
		}
	)
	names(cox_res_atc2) <- atc2names
	
	end.time <- Sys.time()
	end.time
	cat("ELAPSED TIME: ", difftime(end.time, start.time), "\n\n")
	
	plan(sequential)
	
	saveRDS(
		object = cox_res_atc2,
		file = dirty_res_file
	)
} else {
	cox_res_atc2 <- readRDS(dirty_res_file)
}

# tidy and save results - first, for all

tidy_cox_res_atc2_all <- map(
	cox_res_atc2,
	~ broom::tidy(.x$cox_all, conf.int = TRUE)
)
names(tidy_cox_res_atc2_all) <- names(cox_res_atc2)

# remove the empty results
which_empty <- which(map_dbl(tidy_cox_res_atc2_all, nrow) == 0)
if(length(which_empty) != 0){
	tidy_cox_res_atc2_all <- tidy_cox_res_atc2_all[-which_empty]
}

saveRDS(
	object = tidy_cox_res_atc2_all,
	file = paste0(tidy_res_file, "_all.rds")
)

# now for women
tidy_cox_res_atc2_women <- map(
	cox_res_atc2,
	~ broom::tidy(.x$cox_women, conf.int = TRUE)
)
names(tidy_cox_res_atc2_women) <- names(cox_res_atc2)

# remove the empty results
which_empty <- which(map_dbl(tidy_cox_res_atc2_women, nrow) == 0)
if(length(which_empty) != 0){
	tidy_cox_res_atc2_women <- tidy_cox_res_atc2_women[-which_empty]
}

saveRDS(
	object = tidy_cox_res_atc2_women,
	file = paste0(tidy_res_file, "_women.rds")
)

# then for men
tidy_cox_res_atc2_men <- map(
	cox_res_atc2,
	~ broom::tidy(.x$cox_men, conf.int = TRUE)
)
names(tidy_cox_res_atc2_men) <- names(cox_res_atc2)

# remove the empty results
which_empty <- which(map_dbl(tidy_cox_res_atc2_men, nrow) == 0)
if(length(which_empty) != 0){
	tidy_cox_res_atc2_men <- tidy_cox_res_atc2_men[-which_empty]
}

saveRDS(
	object = tidy_cox_res_atc2_men,
	file = paste0(tidy_res_file, "_men.rds")
)
