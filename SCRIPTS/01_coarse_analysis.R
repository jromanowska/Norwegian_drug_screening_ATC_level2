# DESCRIPTION: General analysis of the ATC-level2 results
# AUTHOR: Julia Romanowska
# DATE CREATED: 2022-02-04
# DATE MODIFIED:

# SETUP --------------
library(tidyverse)
library(here)
library(ggrepel)

# READ DATA --------------
# all the results
atc2level_all <- read_csv(
	here(
		"DATA",
		"expose_2_prescr_PD-4prescr-levo-mao-b_age-time-scale_results_atc2level_all.csv"
	)
)
atc2level_all

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

# maybe we don't need _all_ the variables at all times
atc2level_all_compact <- atc2level_all %>%
	select(-(coeff:CI.low), -p.val.group, -users.group, -hover.text)
atc2level_men_compact <- atc2level_men %>%
	select(-(coeff:CI.low), -p.val.group, -users.group, -hover.text)
atc2level_women_compact <- atc2level_women %>%
	select(-(coeff:CI.low), -p.val.group, -users.group, -hover.text)

# ANALYSE --------------
# 1. how many significant per group?
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
	mutate(
		prcnt = ifelse(
			direction_change == "decrease",
			yes = -prcnt,
			no = prcnt
		)
	)
signif_per_group_all

signif_per_group_all %>%
	ggplot(aes(group)) +
	geom_col(aes(x = group, y = prcnt, fill = signif)) +
	geom_label(
		aes(
			x = group, y = 0, label = group
		),
		position = position_dodge(width = 1)
	) +
	# scale_color_manual(values = c("black", "white")) +
	# annotate("label", x = 1, y = 20, label = "non-significant") +
	scale_fill_manual(
		values = c("grey30", "grey70")
	) +
	coord_flip() +
	facet_grid(cols = vars(direction_change), scales = "free_x") +
	theme_minimal() +
	scale_y_continuous(expand = c(0,0)) +
	theme(
		panel.grid = element_blank(),
		axis.title.y = element_blank(),
		axis.text.y = element_blank(),
		panel.spacing.x = unit(0, units = "cm")
	)


