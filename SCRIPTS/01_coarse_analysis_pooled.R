# DESCRIPTION: General analysis of the ATC-level2 results
# AUTHOR: Julia Romanowska
# DATE CREATED: 2022-02-04
# DATE MODIFIED: 2022-08-02

# SETUP --------------
library(tidyverse)
library(here)
library(patchwork)
library(gt)
library(ggiraph)

# READ DATA --------------
# official ATC descriptions + DDD
atc_descr <- read_csv(
	here("DATA", "2020_ATC_Index_with_DDDs_electronic_version.csv")
	) %>%
	rename_with(~ str_replace_all(
		.x,
		pattern = fixed(" "),
		replacement = fixed("_")
	)
	)
atc_descr

# all the results
atc2level_all <- read_csv(
	here(
		"DATA",
		"expose_2_prescr_PD-4prescr-levo-mao-b_age-time-scale_results_atc2level_all.csv"
	)
)
atc2level_all

# maybe we don't need _all_ the variables at all times
atc2level_all_compact <- atc2level_all %>%
	select(-(coeff:CI.low), -p.val.group, -users.group, -hover.text)

# ANALYSE --------------
## 1. how many significant per group? ----
signif_per_group_all <- atc2level_all_compact %>%
	mutate(signif = p.adj.FDR < 0.05) %>%
	mutate(
		direction_change = ifelse(
			HR < 1,
			yes = "decrease",
			no = "increase"
		)
	) %>%
	group_by(group, direction_change) %>%
	count(signif)

drugs_per_group <- atc2level_all_compact %>%
	count(group, name = "drugs_per_group")

signif_per_group_all <- signif_per_group_all %>%
	left_join(drugs_per_group) %>%
	mutate(prcnt = n/drugs_per_group * 100) %>%
	ungroup() %>%
	mutate(group = forcats::fct_rev(as.factor(group))) %>%
	left_join(
		atc_descr %>%
			select(ATC_code, ATC_level_name),
		by = c("group" = "ATC_code")
	) %>%
	mutate(ATC_level_name = stringr::str_to_sentence(ATC_level_name))
signif_per_group_all

# write down for further use
write_delim(
	signif_per_group_all,
	here("DATA", "signif_res_pooled_grouped.txt"),
	delim = "\t"
)

prcnt_signif_incr <- signif_per_group_all %>%
	filter(direction_change == "increase") %>%
	# there were no drugs found to 'increase PD-risk' in L-group!
	tibble::add_row(
		group = "L",
		direction_change = "increase",
		signif = FALSE,
		n = NA,
		drugs_per_group = NA,
		prcnt = 0,
		ATC_level_name = stringr::str_to_sentence(
			atc_descr %>%
				filter(ATC_code == "L") %>%
				pull(ATC_level_name)
		)
	) %>%
	mutate(group = forcats::fct_rev(as.factor(group))) %>%
	ggplot(aes(x = group, y = prcnt, fill = signif)) +
	geom_col(aes()) +
	geom_label(
		aes(
			x = group, y = -5, label = group
		),
		position = position_dodge(width = 1),
		inherit.aes = FALSE
	) +
	scale_fill_manual(
		values = c("grey70", "grey30"),
		name = "significant?",
		labels = c("no", "yes")
	) +
	coord_flip() +
	geom_text(
		aes(label = as.character(n)),
		position = position_stack(vjust = 0.5, reverse = FALSE)
	) +
	theme_minimal() +
	scale_y_continuous(
		expand = c(0,0),
		limits = c(-8, 105),
		breaks = seq(0, 100, 25),
		labels = paste0(seq(0, 100, 25), "%")
	) +
	labs(title = "increasing PD-risk") +
	theme(
		panel.grid = element_blank(),
		axis.title.y = element_blank(),
		axis.title.x = element_blank(),
		axis.text.y = element_blank(),
		panel.spacing.x = unit(0, units = "cm"),
		legend.position = c(0.8, 0.95),
		plot.title = element_text(hjust = 0.3)
	)

prcnt_signif_decr <- signif_per_group_all %>%
	filter(direction_change == "decrease") %>%
	ggplot(aes(x = group, y = -prcnt, fill = signif)) +
	geom_col() +
	scale_fill_manual(
		values = c("grey70", "grey30")
	) +
	coord_flip() +
	geom_text(
		aes(label = as.character(n)),
		position = position_stack(vjust = 0.5, reverse = FALSE)
	) +
	theme_minimal() +
	scale_y_continuous(
		expand = c(0,0),
		limits = c(-105, 0),
		breaks = seq(-100, 0, 25),
		labels = paste0(seq(100, 0, -25), "%")
	) +
	labs(title = "decreasing PD-risk") +
	theme(
		panel.grid = element_blank(),
		axis.title.y = element_blank(),
		axis.title.x = element_blank(),
		axis.text.y = element_blank(),
		panel.spacing.x = unit(0, units = "cm"),
		legend.position = "none",
		plot.margin = unit(c(5.5, 0, 5.5, 5.5), units = "points"),
		plot.title = element_text(hjust = 0.7)
	)

prcnt_signif_decr + prcnt_signif_incr +
	plot_layout(widths = c(0.9, 1))

ggsave(
	here("FIGURES", "prcnt_drugs_per_ATC_group_risk_change.png"),
	width = 9,
	height = 5
)

# which groups are which?
signif_per_group_all <- signif_per_group_all %>%
	left_join(
		atc_descr %>% select(ATC_code, ATC_level_name),
		by = c("group" = "ATC_code")
	)
summary_atc_group_tbl <- signif_per_group_all %>%
		distinct(group, ATC_level_name, drugs_per_group) %>%
		mutate(ATC_level_name = str_to_title(ATC_level_name)) %>%
		select(group, ATC_level_name, drugs_per_group)

knitr::kable(summary_atc_group_tbl)

latex_table <- knitr::kable(
	summary_atc_group_tbl,
	format = "latex"
)

## 2. Ordering of significant results ----
atc2level_signif_compact <- atc2level_all_compact %>%
	filter(p.adj.FDR < 0.05) %>%
	mutate(
		direction_change = ifelse(
			HR < 1,
			yes = "decrease",
			no = "increase"
		),
		name = str_to_title(name)
	) %>%
	arrange(p.adj.FDR)

# save for later use
write_delim(
	atc2level_signif_compact,
	here("DATA", "signif_res_pooled_all.txt"),
	delim = "\t"
)

atc2level_signif_incr <- atc2level_signif_compact %>%
	filter(direction_change == "increase") %>%
	arrange(p.adj.FDR) %>%
	select(ATC_code, name, HR:N.pd)

atc2level_signif_decr <- atc2level_signif_compact %>%
	filter(direction_change == "decrease") %>%
	arrange(p.adj.FDR) %>%
	select(ATC_code, name, HR:N.pd)

# nice latex tables
atc2level_signif_gt <- gt(
		data = atc2level_signif_compact %>%
			arrange(direction_change, HR, ATC_code) %>%
			select(ATC_code, name, N.nonpd, N.pd, HR:p.adj.FDR),
		groupname_col = "direction_change"
	) %>%
	tab_row_group(
		label = "drugs increasing PD risk",
		rows = HR > 1
	) %>%
	tab_row_group(
		label = "drugs decreasing PD risk",
		rows = HR < 1
	) %>%
	cols_merge_range(
		col_begin = HR_l,
		col_end = HR_u
	) %>%
	cols_label(
		ATC_code = "symbol",
		name = "description",
		HR_l = "95% CI",
		p.adj.FDR = "FDR p-value",
		N.nonpd = "non-PD users",
		N.pd = "PD users"
	) %>%
	tab_spanner(
		label = "ATC sub-group",
		columns = c(ATC_code, name)
	) %>%
	tab_spanner(
		label = "Cox hazard risk (HR) estimates",
		columns = HR:p.adj.FDR
	) %>%
	tab_spanner(
		label = "# of users",
		columns = starts_with("N.")
	) %>%
	tab_style(
		style = cell_borders(sides = "right", color = "gray30"),
		locations = cells_body(
			columns = c(name, p.adj.FDR)
		)
	) %>%
	tab_style(
		style = cell_text(style = "italic"),
		locations = cells_body(
			columns = name
		)
	) %>%
	tab_style(
		style = cell_text(weight = "bold"),
		locations = cells_column_labels()
	) %>%
	tab_style(
		style = cell_text(weight = "bold"),
		locations = cells_column_spanners()
	) %>%
	tab_style(
		style = list(
			cell_text(align = "center"),
			cell_fill()
			),
		locations = cells_row_groups()
	) %>%
	fmt_number(
		columns = starts_with("HR")
	) %>%
	fmt_scientific(
		columns = p.adj.FDR
	) %>%
	fmt_number(
		columns = starts_with("N."),
		decimals = 0,
		use_seps = TRUE
	) %>%
	cols_width(
		name ~ px(150)
	) %>%
	tab_footnote(
		footnote = "p-value adjusted for false discovery rate",
		locations = cells_column_labels(columns = p.adj.FDR)
	)
atc2level_signif_gt

atc2level_signif_gt_latex <- as_latex(atc2level_signif_gt)
write_lines(
	atc2level_signif_gt_latex,
	file = here("RESULTS", "signif_res_pooled.tex")
)

## 3. How many per each group per direction of change? ----
gt(
	atc2level_signif_compact %>%
	group_by(direction_change) %>%
	count(group) %>%
	arrange(direction_change, desc(n))
)
# this is not so interesting...

## 4. Plot HR per ATC group ----
ATC_order_by_HR <- atc2level_all_compact %>%
	arrange(HR) %>%
	pull(ATC_code)

plot_all <- ggplot(
	atc2level_all_compact %>%
		mutate(
			ATC_code = 
					 	factor(ATC_code, levels = ATC_order_by_HR, ordered = TRUE),
			hover_text = paste0(
				ATC_code, ", HR: ", HR
			)
		),
	aes(x = HR, y = ATC_code)
	) +
	geom_vline(xintercept = 1) +
	geom_pointrange_interactive(
		aes(xmin = HR_l, xmax = HR_u, color = p.adj.FDR < 0.05, tooltip = hover_text),
		size = 0.5
	) +
	scale_x_continuous(trans = "log") +
	scale_color_manual(values = c("gray80", "gray20")) + # colors for significant vs. non-significant
	theme_minimal()
girafe(ggobj = plot_all, width_svg = 7, height = 9)
