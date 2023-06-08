# DESCRIPTION: Analysis of the sensitivity analyses: 2, 4, and 8 prescriptions
#   required to be cosidered as exposed to a drug
# AUTHOR: Julia Romanowska
# DATE CREATED: 2023-05-22
# DATE MODIFIED: 2023-06-05

# SETUP --------------
library(tidyverse)
library(here)
library(patchwork)
library(gt)

prescr_exposure_strata <- c(2, 4, 8)
results_file_base <- "prescr_PD-4prescr-levo-mao-b_age-time-scale_results_atc2level_time-lag0yrs_all.csv"
prescr_exposure_results_all_files <- paste0(
	"expose_", prescr_exposure_strata, "_", results_file_base
)

names_all_lag_results <- paste0("exposure_", prescr_exposure_strata, "prescr")

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
atc_descr

# all the results
atc2level_all_exposure_prescr <- map(
	prescr_exposure_results_all_files,
	function(cur_file){
	read_csv(
		here("DATA", cur_file)
	)
})
names(atc2level_all_exposure_prescr) <- names_all_lag_results

atc2level_all_exposure_prescr$exposure_2prescr

# maybe we don't need _all_ the variables at all times
atc2level_all_exposure_prescr_compact <- map(names_all_lag_results, function(cur_name){
	cur_results <- atc2level_all_exposure_prescr[[cur_name]]
	cur_results %>%
		select(-(coeff:CI.low), -p.val.group, -users.group, -hover.text) %>%
		add_column(stratum = cur_name)
}) %>%
	bind_rows()

atc2level_all_exposure_prescr_compact <- atc2level_all_exposure_prescr_compact %>%
	mutate(stratum = factor(stratum, levels = names_all_lag_results, ordered = TRUE))

skimr::skim(atc2level_all_exposure_prescr_compact)

atc2level_all_exposure_prescr_compact %>%
	distinct(stratum, ATC_code) %>%
	count(stratum)

orig_signif_ATC_codes <- atc2level_all_exposure_prescr$exposure_2prescr %>%
	filter(p.adj.FDR < 0.05) %>%
	pull(ATC_code)

atc2level_all_exposure_prescr_compact_orig_signif <- atc2level_all_exposure_prescr_compact %>%
	filter(ATC_code %in% orig_signif_ATC_codes)

atc2level_all_exposure_prescr_compact_orig_signif %>%
	distinct(stratum, ATC_code) %>%
	count(stratum)

# DIFFERENCE IN SIGNIFICANCE? ----
atc2level_all_exposure_prescr_compact %>%
	filter(p.adj.FDR < 0.05) %>%
	distinct(ATC_code, name, stratum) %>%
	mutate(
		name = tolower(name),
		ATC_code = paste0(ATC_code, " (", name, ")")
	) %>%
	select(-name) %>%
	add_column(value = TRUE) %>%
	pivot_wider(
		id_cols = "ATC_code",
		names_from = "stratum",
		values_from = "value",
		values_fill = FALSE
	) %>%
	gt(rowname_col = "ATC_code") %>%
	tab_header(
		title = "Comparison of significant findings",
		subtitle = "between various exposure definitions"
	) %>%
	tab_stubhead(label = "ATC code") %>%
	tab_spanner(
		label = "exposure defined by min.",
		columns = names_all_lag_results
	) %>%
	cols_label(
		exposure_2prescr = "2 prescr.",
		exposure_4prescr = "4 prescr."
	) %>%
	cols_width(
		ATC_code ~ px(300)
	) %>%
	tab_style_body(
		style = cell_fill(color = "grey60"),
		fn = function(x) isFALSE(x)
	)

# PLOT ----
clrs <- c("#585CFA", "#63B3F5", "#58FAED")
names(clrs) <- names_all_lag_results

## PLOT ONLY SIGNIFICANT IN ORIGINAL RESULTS ----
## increase risk
ATC_codes_increase <- atc2level_all_exposure_prescr$exposure_2prescr %>%
	filter(p.adj.FDR < 0.05 & HR >= 1) %>%
	pull(ATC_code)
results_4plot_increase <- atc2level_all_exposure_prescr_compact %>%
	filter(ATC_code %in% ATC_codes_increase) %>%
	mutate(
		ATC_code = paste0(ATC_code, " (", tolower(name), ")"),
		stratum = factor(stratum, levels = rev(names_all_lag_results))
	)

## decrease risk
ATC_codes_decrease <- atc2level_all_exposure_prescr$exposure_2prescr %>%
	filter(p.adj.FDR < 0.05 & HR < 1) %>%
	pull(ATC_code)
results_4plot_decrease <- atc2level_all_exposure_prescr_compact %>%
	filter(ATC_code %in% ATC_codes_decrease) %>%
	mutate(
		ATC_code = paste0(ATC_code, " (", tolower(name), ")"),
		stratum = factor(stratum, levels = rev(names_all_lag_results))
	)

plot_prescr_exposure_compare <- function(cur_data, colors = clrs, ncolumns = 5){
	ggplot(cur_data, aes(HR, stratum)) +
		geom_vline(xintercept = 1) +
		geom_linerange(
			aes(xmin = HR_l, xmax = HR_u, colour = stratum),
			linewidth = 1
		) +
		geom_point(
			aes(color = stratum),
			size = 3
		) +
		scale_color_manual(
			values = colors, name = "", breaks = names_all_lag_results
		) +
		scale_x_continuous(trans = "log") +
		facet_grid(
			rows =  vars(ATC_code),
			scales = "free_y",
			switch = "y",
			labeller = label_wrap_gen(50)
		) +
		theme_light() +
		theme(
			axis.title.y = element_blank(),
			axis.text.y = element_blank(),
			axis.ticks.y = element_blank(),
			panel.grid.major.y = element_blank(),
			panel.grid.major.x = element_line(colour = "grey80"),
			strip.text.y.left = element_text(
				angle = 0,
				hjust = 0,
				face = "bold",
				size = 9
			),
			strip.placement = "outside",
			panel.spacing.y = unit(5, units = "points")#,
			#legend.position = "right"
		)
}

plot_signif_incr_risk <- plot_prescr_exposure_compare(cur_data = results_4plot_increase) +
	labs(title = "Drugs originally significantly associated with increased PD risk") +
	theme(legend.position = "bottom")

plot_signif_decr_risk <- plot_prescr_exposure_compare(cur_data = results_4plot_decrease) +
	labs(title = "Drugs originally significantly associated with decreased PD risk") +
	theme(legend.position = "none")

plot_signif_incr_risk / plot_signif_decr_risk +
	plot_layout(
#		guides = "collect",
		heights = c(1, 10/21)
	)
ggsave(
	here("FIGURES", "comparison_sensitivity_results_signif_only.png"),
	height = 12,
	width = 12
)

# ## PLOT ALL ----
# ## increase risk
# ATC_codes_all_increase <- atc2level_all_exposure_prescr$exposure_2prescr %>%
# 	filter(HR >= 1) %>%
# 	pull(ATC_code)
# results_all_4plot_increase <- atc2level_all_exposure_prescr_compact %>%
# 	filter(ATC_code %in% ATC_codes_all_increase) %>%
# 	mutate(
# 		ATC_code = paste0(ATC_code, " (", tolower(name), ")"),
# 		stratum = factor(stratum, levels = rev(names_all_lag_results))
# 	)
# 
# ## decrease risk
# ATC_codes_all_decrease <- atc2level_all_exposure_prescr$exposure_2prescr %>%
# 	filter(HR < 1) %>%
# 	pull(ATC_code)
# results_all_4plot_decrease <- atc2level_all_exposure_prescr_compact %>%
# 	filter(ATC_code %in% ATC_codes_all_decrease) %>%
# 	mutate(
# 		ATC_code = paste0(ATC_code, " (", tolower(name), ")"),
# 		stratum = factor(stratum, levels = rev(names_all_lag_results))
# 	)
# 
# plot_incr_risk <- plot_prescr_exposure_compare(cur_data = results_all_4plot_increase) +
# 	labs(title = "Drugs originally associated with increased PD risk")
# 
# plot_decr_risk <- plot_prescr_exposure_compare(cur_data = results_all_4plot_decrease) +
# 	labs(title = "Drugs originally associated with decreased PD risk")
# 
# plot_incr_risk / plot_decr_risk +
# 	plot_layout(
# 		guides = "collect",
# 		heights = c(1, 7/9)
# 	) & theme(legend.position = "bottom")
# ggsave(
# 	here("FIGURES", "comparison_sensitivity_results_all.png"),
# 	height = 12,
# 	width = 16
# )
