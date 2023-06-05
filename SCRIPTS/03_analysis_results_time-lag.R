# DESCRIPTION: Analysis of the results: time-lagged, ATC level 2, not stratified
# AUTHOR: Julia Romanowska
# DATE CREATED: 2023-05-08
# DATE MODIFIED: 2023-05-28

# SETUP --------------
library(tidyverse)
library(here)
library(patchwork)
library(gt)

time_lag_strata <- c(5, 8, 10)
results_file_base <- "expose_2_prescr_PD-4prescr-levo-mao-b_age-time-scale_results_atc2level"
time_lag_results_all_files <- paste0(
	results_file_base, "_time-lag", time_lag_strata, "yrs_all.csv"
)
orig_results_file <- "expose_2_prescr_PD-4prescr-levo-mao-b_age-time-scale_results_atc2level_time-lag0yrs_all.csv"

names_all_lag_results <- c("orig", paste0("lag", time_lag_strata, "yrs"))

# READ DATA ---------
# official ATC descriptions + DDD
atc_descr <- read_csv(
	here("DATA", "2020_ATC_Index_with_DDDs_electronic_version.csv")
	) %>%
	rename_with(~ str_replace_all(
		.x,
		pattern = fixed(" "),
		replacement = fixed("_")
	)
	) %>%
	mutate(ATC_level_name = stringr::str_to_sentence(ATC_level_name))

# all the results
atc2level_all_time_lags <- map(
	c(orig_results_file, time_lag_results_all_files),
	function(cur_file){
	read_csv(
		here("DATA", cur_file)
	)
})
names(atc2level_all_time_lags) <- names_all_lag_results

atc2level_all_time_lags$orig
atc2level_all_time_lags$lag5yrs

# maybe we don't need _all_ the variables at all times
atc2level_all_time_lags_compact <- map(names_all_lag_results, function(cur_name){
	cur_results <- atc2level_all_time_lags[[cur_name]]
	cur_results %>%
		select(-(coeff:CI.low), -p.val.group, -users.group, -hover.text) %>%
		add_column(stratum = cur_name)
}) %>%
	bind_rows()

atc2level_all_time_lags_compact <- atc2level_all_time_lags_compact %>%
	mutate(stratum = factor(stratum, levels = names_all_lag_results, ordered = TRUE))

skimr::skim(atc2level_all_time_lags_compact)

atc2level_all_time_lags_compact %>%
	distinct(stratum, ATC_code) %>%
	count(stratum)

orig_signif_ATC_codes <- atc2level_all_time_lags$orig %>%
	filter(p.adj.FDR < 0.05) %>%
	pull(ATC_code)

atc2level_all_time_lags_compact_orig_signif <- atc2level_all_time_lags_compact %>%
	filter(ATC_code %in% orig_signif_ATC_codes)

atc2level_all_time_lags_compact_orig_signif %>%
	distinct(stratum, ATC_code) %>%
	count(stratum)

# PLOT ----
clrs <- c("#585CFA", "#63B3F5", "#4EC3DE", "#58FAED")
names(clrs) <- names_all_lag_results

## PLOT ONLY SIGNIFICANT IN ORIGINAL RESULTS ----
## increase risk
ATC_codes_increase <- atc2level_all_time_lags$orig %>%
	filter(p.adj.FDR < 0.05 & HR >= 1) %>%
	pull(ATC_code)
results_4plot_increase <- atc2level_all_time_lags_compact %>%
	filter(ATC_code %in% ATC_codes_increase) %>%
	mutate(ATC_code = paste0(ATC_code, " (", tolower(name), ")"))

## decrease risk
ATC_codes_decrease <- atc2level_all_time_lags$orig %>%
	filter(p.adj.FDR < 0.05 & HR < 1) %>%
	pull(ATC_code)
results_4plot_decrease <- atc2level_all_time_lags_compact %>%
	filter(ATC_code %in% ATC_codes_decrease) %>%
	mutate(ATC_code = paste0(ATC_code, " (", tolower(name), ")"))

plot_time_lag_compare <- function(cur_data, colors = clrs, ncolumns = 5){
	ggplot(cur_data, aes(stratum, HR)) +
		geom_hline(yintercept = 1) +
		geom_linerange(
			aes(ymin = HR_l, ymax = HR_u, colour = stratum),
			linewidth = 1
		) +
		geom_point(
			aes(color = stratum),
			size = 3
		) +
		scale_color_manual(
			values = colors
		) +
		scale_y_continuous(trans = "log") +
		facet_wrap(
			facets = ~ ATC_code,
			ncol = ncolumns,
			scales = "free_y"
		) +
		theme_light() +
		theme(
			axis.title.x = element_blank(),
			axis.text.x = element_blank(),
			strip.text = element_text(hjust = 0)
		)
}

plot_time_lag_compare(cur_data = results_4plot_increase) +
	labs(title = "Drugs originally significantly associated with increased PD risk")
ggsave(
	here("FIGURES", "comparison_time-lag_signif_only_increase_risk.png"),
	width = 12,
	height = 10
)

plot_time_lag_compare(cur_data = results_4plot_decrease) +
	labs(title = "Drugs originally significantly associated with decreased PD risk")
ggsave(
	here("FIGURES", "comparison_time-lag_signif_only_decrease_risk.png"),
	width = 12,
	height = 10
)

### check which are significant when ----
plot_all_signif <- function(cur_data){
	cur_data <- cur_data %>%
		mutate(label = paste0(ATC_code, " (", name, ")"))
	order_by_HR <- cur_data %>%
		arrange(HR) %>%
		pull(label)
	cur_data <- cur_data %>%
		mutate(label = factor(label, levels = order_by_HR, ordered = TRUE))
		
	ggplot(cur_data, aes(HR, label)) +
		geom_vline(xintercept = 1, color = "grey30") +
		geom_pointrange(
			aes(
				xmin = HR_l,
				xmax = HR_u,
				color = group
			)
		) +
		scale_x_continuous(trans = "log") +
		theme_minimal() +
		theme(
			axis.title.y = element_blank()
		)
}

plot_all_signif(atc2level_all_time_lags$orig %>% filter(p.adj.FDR < 0.05)) +
	labs(title = "Time-lag: 0 years")

plot_all_signif(atc2level_all_time_lags$lag5yrs %>% filter(p.adj.FDR < 0.05)) +
	labs(title = "Time-lag: 5 years")

plot_all_signif(atc2level_all_time_lags$lag8yrs %>% filter(p.adj.FDR < 0.05)) +
	labs(title = "Time-lag: 8 years")

plot_all_signif(atc2level_all_time_lags$lag10yrs %>% filter(p.adj.FDR < 0.05)) +
	labs(title = "Time-lag: 10 years")

## PLOT ALL ----
## increase risk
ATC_codes_all_increase <- atc2level_all_time_lags$orig %>%
	filter(HR >= 1) %>%
	pull(ATC_code)
results_all_4plot_increase <- atc2level_all_time_lags_compact %>%
	filter(ATC_code %in% ATC_codes_all_increase) %>%
	mutate(ATC_code = paste0(ATC_code, " (", tolower(name), ")"))

## decrease risk
ATC_codes_all_decrease <- atc2level_all_time_lags$orig %>%
	filter(HR < 1) %>%
	pull(ATC_code)
results_all_4plot_decrease <- atc2level_all_time_lags_compact %>%
	filter(ATC_code %in% ATC_codes_all_decrease) %>%
	mutate(ATC_code = paste0(ATC_code, " (", tolower(name), ")"))

plot_time_lag_compare(cur_data = results_all_4plot_increase) +
	labs(
		title = "Drugs originally associated with increased PD risk"
	)

plot_time_lag_compare(cur_data = results_all_4plot_decrease) +
	labs(
		title = "Drugs originally associated with decreased PD risk"
	)
