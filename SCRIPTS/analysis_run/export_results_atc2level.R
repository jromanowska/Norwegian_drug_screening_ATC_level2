# export results
library(readr)
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)

dir_in <- file.path("..", "PD_screening_ATC_level2")
exposure_strata_names <- as.character(c(2, 4, 8))
time_lag_strata_names <- as.character(c(5, 8, 10))
all_extra_strata_names <- c(time_lag_strata_names, exposure_strata_names)

strata_names <- c("all", "women", "men")

time_lag_results_file_names <- paste0(
	"nice_results_atc_group2_time-lag", time_lag_strata_names, "yrs"
)
exposure_definitions_results_file_names <- paste0(
	"nice_results_atc_group2_exposure", exposure_strata_names, "prescr"
)
all_results_file_base <- c(time_lag_results_file_names,
														exposure_definitions_results_file_names)

# read the ATC codes with their names
all_ATC_codes <- read_delim(
	file.path("..", "..", "Copy_original",
						"ATC_codes", "2020_ATC_Index_with_DDDs_electronic_version.csv"),
	delim = ";", col_names = TRUE
) %>%
	rename(ATC_code = `ATC code`,
				 name = `ATC level name`)

walk(seq_along(all_results_file_base), function(cur_results_number){
	cur_results_file <- all_results_file_base[cur_results_number]
	# set start values
	cur_exposure <- 2
	cur_time_lag <- 0
	
	if(stringr::str_detect(
			string = as.character(cur_results_file),
			pattern = "exposure"
		)){
		cur_exposure <- all_extra_strata_names[cur_results_number]
		cat("cur_exposure: ", cur_exposure, "\n")
	} else {
		cur_time_lag <- all_extra_strata_names[cur_results_number]
		cat("cur_time_lag: ", cur_time_lag, "\n")
	}
	
	# read additional information: number of users who have PD or don't have PD
	all_users_per_atc <- read_delim(
		delim = "\t",
		file.path(
			dir_in, "DATA", paste0(
				"n_users_per_ATC_code2_exposure", cur_exposure,
				"_time-lag", cur_time_lag, "yrs.txt"
			)
		)
	)
	all_users_per_atc <- list(
		women = all_users_per_atc %>%
			filter(sex == 2) %>%
			select(-sex),
		men = all_users_per_atc %>%
			filter(sex == 1) %>%
			select(-sex)
	)
	all_users_per_atc$all <- read_delim(
		delim = "\t",
		file.path(
			dir_in, "DATA",
			paste0(
				"n_users_per_ATC_code2_no_sex_exposure", cur_exposure,
				"_time-lag", cur_time_lag, "yrs.txt"
			)
		)
	)
	
	all_results <- map(strata_names, function(stratum){
		file_in <- paste0(cur_results_file, "_", stratum, ".rds")
		cur_results <- readRDS(file.path(dir_in, "RESULTS", file_in))
		cur_results
		
		cur_users_per_atc <- all_users_per_atc[[stratum]] %>%
			mutate(park_yn = ifelse(park_yn, yes = "N.pd", no = "N.nonpd")) %>%
			pivot_wider(id_cols = ATC_code, names_from = park_yn, values_from = n_users)
		
		# merge
		cur_results <- cur_results %>%
			left_join(all_ATC_codes %>%
									select(ATC_code, name),
								by = "ATC_code") %>%
			mutate(group = substr(ATC_code, 1, 1)) %>%
			mutate(HR = exp(coeff), HR_l = exp(CI.low), HR_u = exp(CI.high)) %>%
			# exclude NO4:
			filter(substr(ATC_code, start = 1, stop = 3) != "N04") %>%
			# exclude unreliable results:
			filter(HR > 0.01) %>%
			# calculate adjusted p-values
			mutate(p.adj.FDR = p.adjust(pvals, method = "BY")) %>%
			# create p-value groups
			mutate(p.val.group = cut(
				p.adj.FDR,
				breaks = c(0, 0.00001, 0.0001, 0.001, 0.01, 0.05, 1),
				labels = c("(0-0.00001]", "(0.00001-0.0001]", "(0.0001-0.001]",
									 "(0.001-0.01]", "(0.01-0.05]", "(0.05-1]"),
				include.lowest = TRUE,
				ordered_result = TRUE)
			) %>%
			left_join(cur_users_per_atc, by = "ATC_code") %>%
			filter(!is.na(N.pd)) %>%
			mutate(users.group = cut(
				N.pd + N.nonpd,
				breaks = c(0, 100, 500, 10000, 30000, 100000, 250000, 1000000, 3000000),
				labels = c("(0-100]", "(100-500]", "(500-10k]", "(10k-30k]", "(30k-100k]",
									 "(100k-250k]", "(250k-1mln]", "(1mln-3mln]"),
				include.lowest = TRUE,
				ordered_result = TRUE)
			) %>%
			# create text for displaying interactively over each point:
			mutate(hover.text = paste0("drug name: ", name,
																 "\nHR (95% CI): ", round(HR, 2),
																 " (", round(HR_l, 2),"; ", round(HR_u, 2),
																 ")\n p-value: ", sprintf("%.1e", pvals),
																 " (adj.: ", sprintf("%.1e", p.adj.FDR), ")",
																 "\n PD-users: ", N.pd,
																 "\n non-PD-users: ", N.nonpd)) %>%
			distinct()
	})
	names(all_results) <- strata_names
	
	readr::write_csv(
		all_results$women,
		paste0(
			"expose_", cur_exposure,
			"_prescr_PD-4prescr-levo-mao-b_age-time-scale_results_atc2level_time-lag",
			cur_time_lag, "yrs_women.csv"
		))
	readr::write_csv(
		all_results$men,
		paste0(
			"expose_", cur_exposure,
			"_prescr_PD-4prescr-levo-mao-b_age-time-scale_results_atc2level_time-lag",
			cur_time_lag, "yrs_men.csv"
		))
	readr::write_csv(
		all_results$all,
		paste0(
			"expose_", cur_exposure,
			"_prescr_PD-4prescr-levo-mao-b_age-time-scale_results_atc2level_time-lag",
			cur_time_lag, "yrs_all.csv"
		))
})

# readr::write_csv(all_results$all, "cur_results_all.csv")
# readr::write_csv(all_results$women, "cur_results_women.csv")
# readr::write_csv(all_results$men, "cur_results_men.csv")
