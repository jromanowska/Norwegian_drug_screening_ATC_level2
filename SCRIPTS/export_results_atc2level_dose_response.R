# export results for the dose-response analysis at ATC level2
library(readr)
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(stringr)

dir_in <- file.path("..", "PD_screening_ATC_level2")
strata_names <- c("all", "women", "men") # but we'll be using only "all"
stratum <- "all"

results_file_base <- paste0(
	"nice_results_atc_group2_time-lag0yrs_dose-response"
)

# read the ATC codes with their names
all_ATC_codes <- read_delim(
	file.path("..", "..", "Copy_original",
						"ATC_codes", "2020_ATC_Index_with_DDDs_electronic_version.csv"),
	delim = ";", col_names = TRUE
) %>%
	rename(ATC_code = `ATC code`,
				 name = `ATC level name`)

# read additional information: person-years and case number per dose
person_years_per_dose <- readRDS(
	file.path(
		dir_in, "DATA", "person_years_and_cases_per_dose_ATC2level.rds"
	)
)
person_years_per_dose$R03

file_in <- paste0(results_file_base, "_", stratum, ".rds")
all_results <- readRDS(file.path(dir_in, "RESULTS", file_in))
all_results$R03

# combine the tables per ATC group
all_tables <- map(
	names(all_results),
	function(atc_group){
		cur_results <- all_results[[atc_group]] %>%
			filter(str_detect(string = term, pattern = "trt2")) %>%
			mutate(
				HR = exp(estimate),
				HR_l = exp(conf.low),
				HR_u = exp(conf.high)
			) %>%
			select(term, HR:HR_u, p.value) %>%
			mutate(
				term = paste0(
					"trt_q", str_sub(term, start = -1)
				)
			) %>%
			add_row( #add the reference row
				term = "trt_q1",
				HR = 1,
				HR_l = NA,
				HR_u = NA,
				p.value = NA
			)
		
		cur_personyears <- person_years_per_dose[[atc_group]]
		personyears_labels <- names(cur_personyears)
		names(cur_personyears) <- c("term", "personyears", "cases")
		
		# merge
		cur_results_table <- left_join(
			cur_personyears,
			cur_results,
			by = "term"
		) %>%
			tibble::add_column(atc_group = atc_group) %>%
			left_join(
				all_ATC_codes %>%
					select(ATC_code, name) %>%
					mutate(name = tolower(name)),
				by = c("atc_group" = "ATC_code")
			)
		return(cur_results_table)
	}
) %>%
	bind_rows()

readr::write_csv(
	all_tables,
	"expose_2_prescr_PD-4prescr-levo-mao-b_age-time-scale_results_atc2level_time-lag0yrs_dose-response_all.csv"
	)

