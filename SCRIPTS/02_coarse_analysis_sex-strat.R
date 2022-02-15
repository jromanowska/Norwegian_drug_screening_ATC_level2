# DESCRIPTION: General analysis of the ATC-level2 results; sex-stratified
# AUTHOR: Julia Romanowska
# DATE CREATED: 2022-02-14
# DATE MODIFIED: 2022-02-14

# SETUP --------------
library(tidyverse)
library(here)
library(patchwork)
library(gt)

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

# sex-stratified
atc2level_men <- read_csv(
	here(
		"DATA",
		"expose_2_prescr_PD-4prescr-levo-mao-b_age-time-scale_results_atc2level_men.csv"
	)
)
atc2level_men
atc2level_women <- read_csv(
	here(
		"DATA",
		"expose_2_prescr_PD-4prescr-levo-mao-b_age-time-scale_results_atc2level_women.csv"
	)
)
atc2level_women

# maybe we don't need _women_ the variables at all times
atc2level_men_compact <- atc2level_men %>%
	select(-(coeff:CI.low), -p.val.group, -users.group, -hover.text)
atc2level_women_compact <- atc2level_women %>%
	select(-(coeff:CI.low), -p.val.group, -users.group, -hover.text)

# ANALYSE --------------
## 1. how many significant per group? ----
### 1A. calculate for women
signif_per_group_women <- atc2level_women_compact %>%
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

drugs_per_group_women <- atc2level_women_compact %>%
	count(group, name = "drugs_per_group")

signif_per_group_women <- signif_per_group_women %>%
	left_join(drugs_per_group_women) %>%
	mutate(prcnt = n/drugs_per_group * 100) %>%
	ungroup() %>%
	mutate(group = as.factor(group))
signif_per_group_women

### 1B. calculate for men
signif_per_group_men <- atc2level_men_compact %>%
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

drugs_per_group_men <- atc2level_men_compact %>%
	count(group, name = "drugs_per_group")

signif_per_group_men <- signif_per_group_men %>%
	left_join(drugs_per_group_men) %>%
	mutate(prcnt = n/drugs_per_group * 100) %>%
	ungroup() %>%
	mutate(group = as.factor(group))
signif_per_group_men

### 1C. join and plot
signif_per_group_all <- bind_rows(
	signif_per_group_women %>%
		add_column(sex = "female"),
	signif_per_group_men %>%
		add_column(sex = "male")
)
signif_per_group_all

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
		sex = c("female", "male")
	) %>%
	# there were no drugs found to 'increase PD-risk' in R and P groups for women!
	tibble::add_row(
		group = c("R", "P"),
		direction_change = "increase",
		signif = FALSE,
		n = NA,
		drugs_per_group = NA,
		prcnt = 0,
		sex = "female"
	) %>%
	ggplot(aes(x = group, y = prcnt, fill = signif)) +
	geom_col(aes()) +
	geom_label(
		aes(x = group, y = -5, label = group),
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
	facet_grid(
		rows = vars(sex)
	) +
	theme_minimal() +
	scale_y_continuous(
		expand = c(0,0),
		limits = c(-8, 105),
		breaks = seq(0, 100, 25),
		labels = paste0(seq(0, 100, 25), "%")
	) +
	scale_x_discrete(
		expand = c(0, 1.2)
	) +
	labs(title = "increasing PD-risk:") +
	theme(
		panel.grid = element_blank(),
		axis.title.y = element_blank(),
		axis.title.x = element_blank(),
		axis.text.y = element_blank(),
		panel.spacing.x = unit(0, units = "cm"),
		legend.position = c(0.8, 0.95),
		plot.title = element_text(hjust = 0.2),
		panel.spacing.y = unit(15, units = "points"),
		strip.text = element_blank()
	)

prcnt_signif_decr <- signif_per_group_all %>%
	filter(direction_change == "decrease") %>%
	# there were no drugs found to 'decrease PD-risk' in B group for women
	#   and in G group for men
	tibble::add_row(
		group = c("B", "G"),
		direction_change = "increase",
		signif = FALSE,
		n = NA,
		drugs_per_group = NA,
		prcnt = 0,
		sex = c("female", "male")
	) %>%
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
	scale_x_discrete(
		expand = c(0, 1.2)
	) +
	facet_grid(
		rows = vars(sex)
	) +
	labs(title = "decreasing PD-risk:") +
	theme(
		panel.grid = element_blank(),
		axis.title.y = element_blank(),
		axis.title.x = element_blank(),
		axis.text.y = element_blank(),
		panel.spacing.x = unit(0, units = "cm"),
		legend.position = "none",
		plot.margin = unit(c(5.5, 0, 5.5, 5.5), units = "points"),
		plot.title = element_text(hjust = 0.8),
		panel.spacing.y = unit(15, units = "points"),
		strip.text = element_blank()
	)

prcnt_signif_decr + prcnt_signif_incr +
	plot_layout(widths = c(0.9, 1)) +
	geom_text(
		data = tibble(
			x = length(levels(signif_per_group_all$group)) + 1,
			y = 0,
			sex = c("female", "male"),
			signif = FALSE
		),
		aes(x = x, y = y, label = paste0(str_to_upper(sex), "S:"))
	) +
	plot_annotation(
		title = "Percentage of drugs in each ATC group that were found\n to change the risk of Parkinson's disease (PD)",
		theme = theme(title = element_text(face = "bold", size = 12))
	)

ggsave(
	here("FIGURES", "prcnt_drugs_per_ATC_group_risk_change_sex_strat.png"),
	width = 9,
	height = 10
)

# 2. Ordering of significant results ----
atc2level_signif_compact <- atc2level_women_compact %>%
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
		select(ATC_code, name, HR:N.pd),
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
# cat(as.character(atc2level_signif_gt_latex))

# 3. How many per each group per direction of change? ----
gt(
	atc2level_signif_compact %>%
	group_by(direction_change) %>%
	count(group) %>%
	arrange(direction_change, desc(n))
)
# this is not so interesting...

# 4. 