# DESCRIPTION: Read the results for level3 in C09: ARBs vs. ACEIs
# AuTHOR: Julia Romanowska
# DATE CREATED: 2023-06-15
# DATE LAST MODIFIED:

# SETUP ----
library(tidyverse)
library(here)

dirty_res_file <- here("RESULTS", "dirty_results_atc_group3_C09.rds")
tidy_res_file <- here("RESULTS", "nice_results_atc_group3_C09")

strata <- c("all", "women", "men")

# READ RESULTS ----
tidy_results_C09 <- map(
	strata,
	function(stratum){
		readRDS(paste0(tidy_res_file, "_", stratum, ".rds")) %>%
			mutate(
				HR = exp(coeff),
				HR_l = exp(CI.low),
				HR_u = exp(CI.high)
			)
	}
)
names(tidy_results_C09) <- strata

tidy_results_C09

walk(
	strata,
	function(stratum){
		cat("\n RESULTS FOR: ", stratum, "\n")
		cur_data <- tidy_results_C09[[stratum]]
		knitr::kable(
			cur_data %>%
				select(group = ATC_code, HR, HR_l, HR_u, pvals)
		) %>% print()
	}
)

