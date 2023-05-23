# DESCRIPTION: Lag-time analyses; setting the end of exposure time to X years
#   before the actual date
#
# AUTHOR: Julia Romanowska
# DATE CREATED: 2023-05-05
# DATE MODIFIED: 2023-05-09

# SETUP ----
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

# check for time lags of 5, 8, and 10
time_lag <- 10
n_parallel_processes <- 3
prescriptions_exposure <- 2

dataset_atc2level_file <- paste0("dataset_ready_for_analysis_exposure",
																 prescriptions_exposure,".rds")

dirty_res_file <- here(
	"RESULTS", paste0("dirty_results_atc_group2_time-lag", time_lag, "yrs.rds")
)
tidy_res_file <- here(
	"RESULTS", paste0("nice_results_atc_group2_time-lag", time_lag, "yrs")
)

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

### PREPARE FOR TIME LAG ----
data_atc2level[
	, time_risk := time_risk - time_lag*365.25
]
# for persons that got their first N04 prescription or died in 2004,
# we need to set the time as 0.5, otherwise tmerge function would complain
data_atc2level <- data_atc2level[time_risk <= 0, time_risk := 0.5 ]

## RUN ANALYSIS -------------
source("analysis.R")
