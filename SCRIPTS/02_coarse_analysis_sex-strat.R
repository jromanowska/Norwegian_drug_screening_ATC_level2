# DESCRIPTION: General analysis of the ATC-level2 results; sex-stratified
# AUTHOR: Julia Romanowska
# DATE CREATED: 2022-02-14
# DATE MODIFIED: 2022-08-16

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

# pooled analysis (to compare afterwards)
atc2level_signif_pooled <- read_delim(
	here("DATA", "signif_res_pooled_all.txt"),
	delim = "\t"
)
atc2level_signif_pooled

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

## 2. Ordering of significant results ----
create_latex_table_signif_res <- function(data_in){
	atc2level_signif_compact <- data_in %>%
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
	
	# nice latex table
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
	
	return(atc2level_signif_gt)
}

atc2level_signif_women_gt_latex <- as_latex(
	create_latex_table_signif_res(atc2level_women_compact)
)
write_lines(
	as.character(atc2level_signif_women_gt_latex),
	file = here("RESULTS", "signif_res_women_nice_table.tex")
)

atc2level_signif_men_gt_latex <- as_latex(
	create_latex_table_signif_res(atc2level_men_compact)
)
write_lines(
	as.character(atc2level_signif_men_gt_latex),
	file = here("RESULTS", "signif_res_men_nice_table.tex")
)


## 3. How many significant findings are significant in both women and men? ----
atc2level_all_compact <- atc2level_men_compact %>%
	tibble::add_column(sex = "male") %>%
	bind_rows(
		atc2level_women_compact %>%
		tibble::add_column(sex = "female")
	)
	
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
atc2level_signif_compact

duplicated_res <- atc2level_signif_compact %>%
	count(ATC_code) %>%
	filter(n > 1) %>%
	distinct(ATC_code) %>%
	pull()
duplicated_res
length(duplicated_res)

knitr::kable(
	atc2level_signif_compact %>%
		filter(ATC_code %in% duplicated_res) %>%
		arrange(ATC_code, sex, direction_change)
)

# check whether any of these results change direction
atc2level_signif_compact %>%
	filter(ATC_code %in% duplicated_res) %>%
	arrange(ATC_code, sex, direction_change) %>%
	group_by(ATC_code) %>%
	count(direction_change)
# NOPE!

knitr::kable(
	atc2level_signif_compact %>%
		filter(ATC_code %in% duplicated_res) %>%
		distinct(ATC_code, direction_change) %>%
		count(direction_change),
	caption = "How many significant in each direction?"
)

# which where decreasing, which increasing?
atc2level_signif_compact %>%
	filter(ATC_code %in% duplicated_res) %>%
	arrange(direction_change, ATC_code) %>%
	distinct(ATC_code, direction_change)

## 4. How the significant results from pooled analysis compare to the sex-strat? ----
atc2level_signif_compare_pool_strat <- atc2level_signif_pooled %>%
	left_join(
		atc2level_all_compact %>%
			select(ATC_code, pvals, HR:p.adj.FDR, sex) %>%
			pivot_wider(
				id_cols = ATC_code,
				names_from = sex,
				values_from = pvals:p.adj.FDR
			),
		by = "ATC_code"
	)

# save for late use
write_delim(
	atc2level_signif_compare_pool_strat,
	file = here("DATA", "signif_res_compare_pool_strat.txt"),
	delim = "\t"
)

strat_compare_pool_gt <- gt(
	data = atc2level_signif_compare_pool_strat %>% 
	  # arrange according to pooled results:
		arrange(direction_change, HR) %>%
		select(-starts_with("pvals"), -(group:N.pd)) %>%
		# reorder columns:
		select(ATC_code, name, direction_change,
					 ends_with("_male"), ends_with("_female")),
	groupname_col = "direction_change"
	) %>%
	cols_merge_range(
		col_begin = HR_l_male,
		col_end = HR_u_male
	) %>%
	cols_merge_range(
		col_begin = HR_l_female,
		col_end = HR_u_female
	) %>%
	cols_label(
		ATC_code = "symbol",
		name = "description",
		HR_male = "HR",
		HR_female = "HR",
		HR_l_male = "95% CI",
		HR_l_female = "95% CI",
		p.adj.FDR_male = "FDR p-value",
		p.adj.FDR_female = "FDR p-value"
	) %>%
	tab_spanner(
		label = "males",
		columns = HR_male:p.adj.FDR_male
	) %>%
	tab_spanner(
		label = "females",
		columns = HR_female:p.adj.FDR_female
	)%>%
	tab_spanner(
		label = "ATC sub-group",
		columns = c(ATC_code, name)
	) %>%
	tab_spanner(
		label = "Cox hazard risk (HR) estimates",
		columns = HR_male:p.adj.FDR_female
	)  %>%
	tab_style(
		style = cell_borders(sides = "right", color = "gray30"),
		locations = cells_body(
			columns = c(name, p.adj.FDR_male)
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
		columns = starts_with("p.adj.FDR")
	) %>%
	cols_width(
		name ~ px(150)
	) %>%
	tab_footnote(
		footnote = "p-value adjusted for false discovery rate",
		locations = cells_column_labels(columns = starts_with("p.adj.FDR"))
	)
strat_compare_pool_gt

strat_compare_pool_latex <- as_latex(
	strat_compare_pool_gt
)
write_lines(
	as.character(strat_compare_pool_latex),
	file = here("RESULTS", "signif_res_compare_sex_nice_table.tex")
)

### 4A. Plot the comparison ----
atc2level_signif_compare_pool_strat_4plotting <-
	atc2level_signif_compare_pool_strat %>% 
	  # arrange according to pooled results:
		arrange(direction_change, HR) %>%
		# reorder columns:
		select(ATC_code, name, direction_change,
					 -starts_with("pvals"), HR, HR_l, HR_u, p.adj.FDR,
					 -N.nonpd, -N.pd,
					 ends_with("_male"), ends_with("_female"),
					 group
					 ) %>%
		mutate(ATC_code = factor(
			ATC_code, levels = sort(unique(ATC_code), decreasing = TRUE)
		))
atc2level_signif_compare_pool_strat_4plotting

nice_colors <- sanzo::sanzo.trio("c232")
names_categ <- c("not stratified", "females only", "males only")
names(nice_colors) <- names_categ

ggplot(
	atc2level_signif_compare_pool_strat_4plotting %>%
		pivot_longer(
			cols = c(HR, HR_male, HR_female),
			names_to = "estimate_cat",
			values_to = "HR"
		),
	aes(HR, ATC_code)
) +
	geom_vline(xintercept = 1, color = "gray40", lwd = 0.8) +
	geom_segment(
		data = atc2level_signif_compare_pool_strat_4plotting,
		aes(
			x = HR_male,
			xend = HR_female,
			y = ATC_code,
			yend = ATC_code
		),
		color = as.character(nice_colors[1])
	) +
	geom_point(aes(color = estimate_cat), size = 1.8) +
	scale_color_manual(
		values = as.character(nice_colors),
		name = "analysis type",
		labels = names_categ
	) +
	scale_x_continuous(
		trans = "log",
		breaks = c(0.3, 0.6, 1.0, 1.5, 2.0),
		labels = as.character(c(0.3, 0.6, 1.0, 1.5, 2.0))
	) +
	xlab("HR") +
	ylab("ATC group") +
	theme_minimal()

ggsave(
	here("FIGURES", "compare_estimates_pooled_sex-strat_signif.png")
)

### 4B. Plot for Julia's poster ----
ATC_order_by_HR <- atc2level_signif_pooled %>%
	arrange(HR) %>%
	pull(ATC_code)

label <- atc2level_signif_pooled %>%
	distinct(ATC_code, name) %>%
	mutate(label = paste0(ATC_code, " (", name, ")")) %>%
	pull(label)
names(label) <- atc2level_signif_pooled %>%
	distinct(ATC_code) %>%
	pull(ATC_code)

atc2level_signif_compare_pool_strat_4poster <- atc2level_signif_pooled %>%
	tibble::add_column(stratum = "pooled") %>%
	bind_rows(
		atc2level_all_compact %>%
			filter(ATC_code %in% ATC_order_by_HR) %>%
			rename(stratum = sex)
	) %>%
	mutate(
		signif = p.adj.FDR < 0.05,
		ATC_code = factor(ATC_code, levels = ATC_order_by_HR, ordered = TRUE)
	)
atc2level_signif_compare_pool_strat_4poster

compare_colors <- c("salmon", "navy", "gray30")

ggplot(
	atc2level_signif_compare_pool_strat_4poster,
	aes(HR, stratum)
	) +
	geom_vline(xintercept = 1) +
	geom_pointrange(aes(xmin = HR_l, xmax = HR_u, color = stratum), size = 1.3) +
	facet_grid(
		rows = vars(ATC_code),
		labeller = as_labeller(label)
		#labeller = labeller(label = label_value(vars(label)))
	) +
	scale_x_continuous(
		trans = "log",
		breaks = c(0.2, 0.4, 0.8, 1, 1.3, 1.6, 2.1)
	) +
	scale_color_manual(values = compare_colors) +
	theme_minimal() +
	theme(
		axis.title.y = element_blank(),
		axis.text.y = element_blank(),
		strip.text.y = element_text(angle = 0, hjust = 0),
		legend.position = c(0.2, 0.2),
		text = element_text(size = 30)
	)
ggsave(
	filename = here("FIGURES", "pooled_compare_sex_4poster.png"),
	# filename = here("FIGURES", "pooled_compare_sex_4poster.svg"),
	height = 20, width = 25
)
