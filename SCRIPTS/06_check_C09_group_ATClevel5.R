# TRYING META-ANALYSIS ON C09 GROUP
# AUTHOR: Julia Romanowska
# DATE CREATED: 2023-06-09
# DATE LAST MODIFIED: 

# SETUP ----
library(tidyverse)
library(here)
library(metafor)

ARB_codes <- c("C09A", "C09B")
ACE_codes <- c("C09C", "C09D")

# READ DATA ----
(
	atc_level5_results <- read_csv(
		here("..", "RESULTS_RAW", "expose_2_prescr_PD-4prescr-levo-mao-b_age-time-scale_results_Magne.csv")
	)
)

# FILTER ----
(c09_ARBs <- atc_level5_results %>%
	filter(str_starts(ATC_code, paste(ARB_codes, collapse = "|")))
)

(c09_ACEs <- atc_level5_results %>%
	filter(str_starts(ATC_code, paste(ACE_codes, collapse = "|")))
)

# ANALYSIS ----
(meta_HR_ACEs <- rma(koeff, sei = SE, data = c09_ACEs, method = "REML"))

(meta_HR_ARBs <- rma(koeff, sei = SE, data = c09_ARBs, method = "REML"))

exp(as.numeric(meta_HR_ACEs$beta))
exp(as.numeric(meta_HR_ACEs$se))

exp(as.numeric(meta_HR_ARBs$beta))
exp(as.numeric(meta_HR_ARBs$se))

