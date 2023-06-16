# DESCRIPTION: General analysis of the ATC-level2 results
# AUTHOR: Julia Romanowska
# DATE CREATED: 2023-06-07
# DATE MODIFIED: 2023-06-16

# SETUP --------------
library(tidyverse)
library(here)
library(gt)
library(flextable)
library(data.table)
library(patchwork)

dose_levels <- c("no use", ">2 prescr.", "25th percentile", "50th percentile",
								 "75th percentile", "90th percentile")

# READ DATA --------------
dose_response_file <- "expose_2_prescr_PD-4prescr-levo-mao-b_age-time-scale_results_atc2level_time-lag0yrs_dose-response_all.csv"
(dose_response <- read_csv(
	here("DATA", dose_response_file)
))

(atc2level_signif_compact <- read_delim(
	here("DATA", "signif_res_pooled_all.txt"),
	delim = "\t"
))
signif_results_atc <- atc2level_signif_compact$ATC_code

# CREATE NICE TABLES ----
## SIGNIFICANT RESULTS ONLY ----
(
	dose_response_signif_only_for_gt <- dose_response %>%
	filter(atc_group %in% signif_results_atc) %>%
	mutate(
		atc_group = paste0(
			atc_group, " (", name, ")"
		),
		estimates = if_else(
			HR == 1,
			"1.0 (ref)",
			paste0(
			signif(HR, 2), " (", signif(HR_l, 2), "-", signif(HR_u, 2), ")"
			)
		),
		p.value = signif(p.value, 2)
	) %>%
	select(atc_group, term, personyears, cases, estimates, p.value)
)

dose_response_signif_only_for_gt <- dose_response_signif_only_for_gt %>%
  	mutate(
  		term = case_when(
  			term == "trt_q1" ~ dose_levels[1],
  			term == "trt_q2" ~ dose_levels[2],
  			term == "trt_q3" ~ dose_levels[3],
  			term == "trt_q4" ~ dose_levels[4],
  			term == "trt_q5" ~ dose_levels[5],
  			term == "trt_q6" ~ dose_levels[6]
  		)
  	)

(
dose_response_table_gt <- gt(
	dose_response_signif_only_for_gt,
	groupname_col = "atc_group",
	rowname_col = "term"
) %>%
		fmt_number(
			columns = c("personyears", "cases"), decimals = 0
		) %>%
		fmt_scientific(
			columns = "p.value"
		)
)

writeLines(
	as_latex(dose_response_table_gt),
	here("RESULTS", "signif_res_dose_response.tex")
)

# nicer - with flextable
small_border <- officer::fp_border(color = "gray", width = 1)
header_border <- officer::fp_border(color = "gray30", width = 2)

dose_response_table_flex <- as_grouped_data(
  dose_response_signif_only_for_gt %>%
  	mutate(
  		p.value = if_else(
  			is.na(p.value),
  			"",
  			sprintf("%.2e", p.value)
  		)
  	),
  groups = "atc_group"
) %>%
  as_flextable() %>% 
  bold(bold = TRUE, part = "header") %>%
  align(i = ~ !is.na(atc_group), align = "center") %>% 
  bold(i = ~ !is.na(atc_group)) %>%
  set_header_labels(
    estimates = "HR (95% CI)",
    term = "dose stratum"
  ) %>%
  border_remove() %>%
  hline_bottom(
    j = 1:5,
    border = header_border
  ) %>%
  bg(
    i = ~!is.na(atc_group),
    bg = "gray80"
  ) %>%
  hline_top(
    j = 1:5,
    border = header_border
  ) %>%
  border_inner_h(
    border = small_border,
    part = "body"
  ) %>%
  fix_border_issues() %>%
  width(
    j = "estimates",
    width = 5,
    unit = "cm"
  )
dose_response_table_flex %>%
  width(
    j = "p.value",
    width = 3,
    unit = "cm"
  )


# FOREST PLOT ----
## function to plot nice estimates + CI ----
hr_plot <- function(cur_data, x_range){
  cur_data %>%
    filter(!is.na(HR_l)) %>%
    ggplot(aes(HR, 1)) +
    geom_vline(xintercept = 1, color = "grey30") +
    geom_pointrange(aes(xmin = HR_l, xmax = HR_u)) +
    #scale_x_continuous(breaks = c(0, 1, 2.5, 5, 7.5, 10)) +
    coord_cartesian(expand = FALSE, xlim = x_range) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank()
    )
}

## create table with plots ----
tmp_data_table <- as.data.table(
  dose_response %>% 
    filter(atc_group %in% signif_results_atc) %>%
    select(atc_group, term, starts_with("HR")) %>%
  	filter(!is.na(HR_l)) %>%
    arrange(atc_group, term)
)

x_range <- range(tmp_data_table$HR)
x_range

plots_per_atc_group <- tmp_data_table[
  , list(gg_plot = list(hr_plot(.SD, x_range))), by = c("atc_group", "term")
]

# where to add the plots
n_atc_group <- length(unique(plots_per_atc_group$atc_group))
n_terms <- length(unique(plots_per_atc_group$term))
idx_atc_rows <- seq.int(from = 1, by = n_terms + 1, length.out = n_atc_group)
idx_add_plot <- map(
  idx_atc_rows,
  function(cur_atc_row){
    seq.int(from = cur_atc_row + 2, by = 1, length.out = n_terms - 1)
  }
) %>% do.call(what = c, args = .)

(
	dose_resp_table_with_plots <- dose_response_table_flex %>%
  width(
    j = "p.value",
    width = 7,
    unit = "cm"
  ) %>%
  append_chunks(
    i = idx_add_plot,
    j = "p.value",
    gg_chunk(
      value = plots_per_atc_group$gg_plot,
      width = 2,
      height = 0.8
    )
  )
)
save_as_html(
	dose_resp_table_with_plots,
	path = here("RESULTS", "dose_response_table_with_plots.html")
)

# ONLY PLOT OF CHOSEN RESULTS ----
dose_response_signif_only <- dose_response %>%
	filter(!is.na(HR_l)) %>%
	filter(term != "trt_q2") %>%
	mutate(
		term = case_when(
  			term == "trt_q3" ~ dose_levels[3],
  			term == "trt_q4" ~ dose_levels[4],
  			term == "trt_q5" ~ dose_levels[5],
  			term == "trt_q6" ~ dose_levels[6]
  		),
		ATC_code = atc_group,
		atc_group = paste0(
			atc_group, " (", name, ")"
		)
	) %>%
	select(ATC_code, atc_group, term, starts_with("HR"), p.value) %>%
	mutate(
		term = factor(
			term, levels = dose_levels[-1], ordered = TRUE
		)
	)


clrs <- c("#4E77DE", "#63B3F5", "#4EC3DE", "#58FAED")
names(clrs) <- dose_levels[-(1:2)]

## increase risk
ATC_codes_increase <- atc2level_signif_compact %>%
	filter(HR > 1.00) %>%
	pull(ATC_code)
results_4plot_increase <- dose_response_signif_only %>%
	filter(ATC_code %in% ATC_codes_increase)

## decrease risk
ATC_codes_decrease <- atc2level_signif_compact %>%
	filter(HR < 1.00) %>%
	pull(ATC_code)
results_4plot_decrease <- dose_response_signif_only %>%
	filter(ATC_code %in% ATC_codes_decrease)

plot_dose_resp_compare <- function(cur_data, colors = clrs, ncolumns = 5){
	ggplot(cur_data, aes(term, HR)) +
		geom_hline(yintercept = 1) +
		geom_linerange(
			aes(ymin = HR_l, ymax = HR_u, colour = term),
			linewidth = 1
		) +
		geom_point(
			aes(color = term),
			size = 2
		) +
		scale_color_manual(
			values = colors,
			name = "# of prescriptions"
		) +
		scale_y_continuous(
			trans = "log",
			labels = function(x){round(x, 2)}
		) +
		# coord_cartesian(
		# 	ylim = ylim
		# ) +
		facet_wrap(
			facets = ~ atc_group,
			ncol = ncolumns,
			# scales = "free_y",
			labeller = label_wrap_gen(30)
		) +
		theme_light() +
		theme(
			axis.title.x = element_blank(),
			axis.text.x = element_blank(),
			axis.ticks.x = element_blank(),
			panel.grid.major.x = element_blank(),
			panel.grid.major.y = element_line(colour = "grey80"),
			panel.spacing.y = unit(5, units = "points"),
			strip.text = element_text(hjust = 0, face = "bold")#,
			#legend.position = "right"
		)
}

plot_signif_incr_risk <- plot_dose_resp_compare(
		cur_data = results_4plot_increase
	) +
	labs(title = "Drugs originally significantly associated with increased PD risk") +
	theme(legend.position = "bottom", legend.justification = c(0,0))

plot_signif_decr_risk <- plot_dose_resp_compare(
		cur_data = results_4plot_decrease
		) +
	labs(title = "Drugs originally significantly associated with decreased PD risk") +
	theme(legend.position = "none")

plot_signif_decr_risk / plot_signif_incr_risk +
	plot_layout(
		heights = c(2/5, 1)
	)
ggsave(
	here("FIGURES", "comparison_dose_response_results_signif_only.png"),
	height = 11,
	width = 12
)

