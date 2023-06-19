# DESCRIPTION: Testing interaction between drug usage and sex: R05, G04 and G03
#
# AUTHOR: Julia Romanowska
# DATE CREATED: 2023-02-13
# DATE MODIFIED:

# SETUP ----------------------

library(data.table)
library(tidyverse)
library(lubridate)
library(here)
library(broom)
library(survival)

dataset_atc2level_file <- "dataset_ready_for_analysis.rds"
dirty_res_file <- here("RESULTS", "dirty_results_sex_interaction_atc_group2.rds")
tidy_res_file <- here("RESULTS", "nice_results_sex_interaction_atc_group2")

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
atc2names_all <- colnames(data_atc2level[1:2, A01:V08])
atc2names_all

other_colnames <- setdiff(
	colnames(data_atc2level),
	atc2names_all
)
other_colnames

atc2names_chosen <- c("G03", "G04", "N07", "R05")

cox_res_sex_interact_atc2 <- map(
	atc2names_chosen,
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
		
		cur_data <- tmerge(
			data1 = cur_data,
			data2 = cur_data,
			id = indiv_ID,
			park = event(time_risk, park_yn),
			trt = tdc(mtid)
		)
		
		cur_data <- cur_data %>% 
			mutate(tstart = if_else(is.na(mtid), 
															tstart + age0, tstart))
		
		cur_data <- cur_data %>% 
			mutate(tstart = if_else(tstart == 0, 
															tstart + age0, tstart))
		
		cur_cox_res_all <- coxph(
			Surv(tstart, tstop, park) ~ 
				as.factor(edu_level) + sex*trt,
			data = cur_data
		)
		
		return(cur_cox_res_all)
	}
)
names(cox_res_sex_interact_atc2) <- atc2names_chosen

saveRDS(cox_res_sex_interact_atc2, dirty_res_file)

cox_res_sex_interact_atc2
