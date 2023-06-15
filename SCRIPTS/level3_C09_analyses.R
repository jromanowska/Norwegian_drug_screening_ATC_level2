# DESCRIPTION: Cox regression on per ATC-group data
#    S:/Project/DRONE/Magne/test4/alder_kvinner.R
#
# AUTHOR: Julia Romanowska
# DATE CREATED: 2023-06-14
# DATE MODIFIED: 2023-06-15

# SETUP ----------------------

library(data.table)
library(stringr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(here)
library(purrr)
library(future)
library(future.apply)
library(broom)
library(survival)

n_parallel_processes <- 2
prescriptions_exposure <- 2
time_lag <- 0

dataset_atc2level_file <- paste0("dataset_ready_for_analysis_exposure",
																 prescriptions_exposure,
																 "_level3_C09.rds")

dirty_res_file <- here("RESULTS", "dirty_results_atc_group3_C09.rds")
tidy_res_file <- here("RESULTS", "nice_results_atc_group3_C09")

## READ DATA -----------------
data_atc2level <- readRDS(
	here("DATA", dataset_atc2level_file)
)
data_atc2level

dim(data_atc2level)
data_atc2level[,.N, by = park_yn]

## FILTER -----------
# remove all persons who were 23 or younger in 2004
data_atc2level <- data_atc2level[age >= 24]
dim(data_atc2level)
data_atc2level[,.N, by = park_yn]

# remove all that got at least one prescription of N04 in 2004
#   THIS WAS DONE IN PREPARATION!

# categorise age into 5-year intervals
data_atc2level[,
							 age_group := cut(age, breaks = c(seq(25, 85, by = 5), 115), 
							 								 right = FALSE)
]
data_atc2level

# calculate age (in days!) at the time of start of follow up
data_atc2level[,
							 age0 := (age * 365.25) + 365.25/2
]

## RUN ANALYSIS -------------
# get all the names of drug groups
atc2names <- c("ACEI", "ARB")

other_colnames <- setdiff(
	colnames(data_atc2level),
	atc2names
)
other_colnames

if(!file.exists(dirty_res_file)){
	# prepare for parallelization
	options(future.globals.maxSize = 3 *10^9 * 1024^2)
	plan(multisession, workers = n_parallel_processes)
	
	start.time <- Sys.time()
	start.time
	
	cox_res_atc2 <- future_lapply(
		atc2names,
		function(atc){
			cur_columns <- c(other_colnames, atc)
			cur_data <- as_tibble(
				data_atc2level
			) %>%
				select(all_of(cur_columns)) %>%
				rename(mtid = all_of(atc)) %>%
				# only those who had the first prescription _after_ 2004
				filter(time_risk > 0) %>%
				# fix problems with negative time
				mutate(mtid = if_else(mtid < 0, 0.5, mtid)) %>%
				# age at exposure
				mutate(mtid = mtid + age0) %>%
				# age at PD onset or censoring event
				mutate(time_risk = time_risk + age0)
			
			if(time_lag != 0){
				cur_data <- cur_data %>%
					mutate(time_risk = time_risk - time_lag*365.25)
			}
			cur_data <- tmerge(
				data1 = cur_data,
				data2 = cur_data,
				id = indiv_ID,
				park = event(time_risk, park_yn),
				trt = tdc(mtid)
			)
			
			cur_data <- cur_data %>% 
				mutate(tstart = if_else(is.na(mtid), 
																tstart + age0, tstart)) %>% 
				mutate(tstart = if_else(tstart == 0, 
																tstart + age0, tstart))
			
			cur_cox_res_all <- coxph(
				Surv(tstart, tstop, park) ~ 
					trt + as.factor(edu_level) + as.factor(sex),
				data = cur_data
			)
			
			cur_cox_res_women <- coxph(
				Surv(tstart, tstop, park) ~ 
					trt + as.factor(edu_level),
				data = cur_data %>% filter(sex == 2)
			)
			
			cur_cox_res_men <- coxph(
				Surv(tstart, tstop, park) ~ 
					trt + as.factor(edu_level),
				data = cur_data %>% filter(sex == 1)
			)
			
			return(
				list(
					all = cur_cox_res_all,
					women = cur_cox_res_women,
					men = cur_cox_res_men
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
	~ broom::tidy(.x$all, conf.int = TRUE)
)
names(tidy_cox_res_atc2_all) <- names(cox_res_atc2)

# remove the empty results
which_empty <- which(map_dbl(tidy_cox_res_atc2_all, nrow) == 0)
if(length(which_empty) != 0){
	tidy_cox_res_atc2_all <- tidy_cox_res_atc2_all[-which_empty]
}

nice_results <-
	tibble(
		ATC_code = names(tidy_cox_res_atc2_all),
		pvals = map_dbl(tidy_cox_res_atc2_all, pluck, "p.value", 1),
		coeff = map_dbl(tidy_cox_res_atc2_all, pluck, "estimate", 1),
		z.obs = map_dbl(tidy_cox_res_atc2_all, pluck, "statistic", 1),
		SE = map_dbl(tidy_cox_res_atc2_all, pluck, "std.error", 1),
		CI.high = map_dbl(tidy_cox_res_atc2_all, pluck, "conf.high", 1),
		CI.low = map_dbl(tidy_cox_res_atc2_all, pluck, "conf.low", 1)
	)

saveRDS(
	object = nice_results,
	file = paste0(tidy_res_file, "_all.rds")
)

# now for women
tidy_cox_res_atc2_women <- map(
	cox_res_atc2,
	~ broom::tidy(.x$women, conf.int = TRUE)
)
names(tidy_cox_res_atc2_women) <- names(cox_res_atc2)

# remove the empty results
which_empty <- which(map_dbl(tidy_cox_res_atc2_women, nrow) == 0)
if(length(which_empty) != 0){
	tidy_cox_res_atc2_women <- tidy_cox_res_atc2_women[-which_empty]
}

nice_results <-
	tibble(
		ATC_code = names(tidy_cox_res_atc2_women),
		pvals = map_dbl(tidy_cox_res_atc2_women, pluck, "p.value", 1),
		coeff = map_dbl(tidy_cox_res_atc2_women, pluck, "estimate", 1),
		z.obs = map_dbl(tidy_cox_res_atc2_women, pluck, "statistic", 1),
		SE = map_dbl(tidy_cox_res_atc2_women, pluck, "std.error", 1),
		CI.high = map_dbl(tidy_cox_res_atc2_women, pluck, "conf.high", 1),
		CI.low = map_dbl(tidy_cox_res_atc2_women, pluck, "conf.low", 1)
	)

saveRDS(
	object = nice_results,
	file = paste0(tidy_res_file, "_women.rds")
)

# then for men
tidy_cox_res_atc2_men <- map(
	cox_res_atc2,
	~ broom::tidy(.x$men, conf.int = TRUE)
)
names(tidy_cox_res_atc2_men) <- names(cox_res_atc2)

# remove the empty results
which_empty <- which(map_dbl(tidy_cox_res_atc2_men, nrow) == 0)
if(length(which_empty) != 0){
	tidy_cox_res_atc2_men <- tidy_cox_res_atc2_men[-which_empty]
}

nice_results <-
	tibble(
		ATC_code = names(tidy_cox_res_atc2_men),
		pvals = map_dbl(tidy_cox_res_atc2_men, pluck, "p.value", 1),
		coeff = map_dbl(tidy_cox_res_atc2_men, pluck, "estimate", 1),
		z.obs = map_dbl(tidy_cox_res_atc2_men, pluck, "statistic", 1),
		SE = map_dbl(tidy_cox_res_atc2_men, pluck, "std.error", 1),
		CI.high = map_dbl(tidy_cox_res_atc2_men, pluck, "conf.high", 1),
		CI.low = map_dbl(tidy_cox_res_atc2_men, pluck, "conf.low", 1)
	)

saveRDS(
	object = nice_results,
	file = paste0(tidy_res_file, "_men.rds")
)
