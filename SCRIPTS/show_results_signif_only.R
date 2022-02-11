# Create a pretty visualization of the results

# 1. SETUP ----
library(tidyverse)
library(patchwork)
library(here)

# 2. READ DATA ----

#  2A. COX RESULTS ----

# n_users_cox <- 

cox_results <- read_csv(
		here("..", "RESULTS_RAW",
			"expose_2_prescr_PD-4prescr-levo-mao-b_age-time-scale_results_Magne.csv"
		)
	)
cox_results

p.adj.thresh <- 0.05
# merge and filter
cox_results_signif <- cox_results %>%
	# exclude NO4:
	filter(substr(ATC_code, start = 1, stop = 3) != "N04") %>%
	# exclude unreliable results:
	filter(HR > 0.01) %>%
	# calculate adjusted p-values
	filter(!is.na(N.pd)) %>%
	# filter the most significant results decreasing the PD risk:
	filter(HR < 1 & p.adj.FDR <= p.adj.thresh) %>%
	distinct()

	# # create a nice label
	# mutate(label_descr = sprintf(
	# 	"%s (%s), adj.p-val: %.1e, no.PD-users: %i, no.non-PD-users: %i",
	# 	ATC_code, name, p.adj.FDR, N.pd, N.nonpd
	# ))

cox_results_signif

#  2B. TIME-LAG RESULTS ----



# 3. ANALYSE ----

# how many in each ATC group?
stats_atc_group <- cox_results_signif %>%
	count(group) %>%
	arrange(desc(n))

knitr::kable(stats_atc_group)


# 4. PLOT ----

# order the results by p-value:
order_pval <- sort(
	cox_results_signif$p.adj.FDR,
	decreasing = TRUE,
	index.return = TRUE)$ix
cox_results_signif$ATC_code <- factor(
	x = cox_results_signif$ATC_code,
	levels = cox_results_signif$ATC_code[order_pval],
	ordered = TRUE
)

y_labels <- map_chr(
	range(cox_results_signif$p.adj.FDR),
	function(x){
		sprintf("%.1e", x)
	})

plot_HR <- ggplot(cox_results_signif, aes(HR, ATC_code)) +
	geom_errorbarh(aes(xmin = HR_l, xmax = HR_u, col = group)) +
	geom_point(aes(col = group), size = 2) +
	scale_x_continuous(trans = "log") +
	# showing only 1st and last value on y-axis
	scale_y_discrete(
		breaks = levels(cox_results_signif$ATC_code)[c(
				1,
				length(levels(cox_results_signif$ATC_code))
			)],
		labels = rev(y_labels)
	) +
	scale_colour_discrete("ATC group") +
	xlab("log(Hazard Ratio)") +
	ylab("FDR-adj. p-value") +
	theme(
		axis.text.y = element_text(hjust = 0),
		plot.title.position = "plot",
		plot.title = element_text(hjust = 0.5),
		plot.subtitle = element_text(hjust = 0.5, face = "italic")
	)

plot_n_PD <- ggplot(cox_results_signif, aes(N.pd/1000, ATC_code)) +
	geom_point(aes(col = group), size = 2) +
	scale_color_discrete(guide = "none") +
	xlab("# PD-users\n (thousands)") +
	theme(
		axis.text.y = element_blank(),
		axis.title.y = element_blank()
	)

plot_n_non_PD <- ggplot(cox_results_signif, aes(N.nonpd/1e06, ATC_code)) +
	geom_point(aes(col = group), size = 2) +
	scale_color_discrete(guide = "none") +
	xlab("# non-PD-users\n (millions)") +
	theme(
		axis.text.y = element_blank(),
		axis.title.y = element_blank()
	)

plot_HR + plot_n_PD + plot_n_non_PD +
	plot_annotation(
		title = "Drugs that significantly reduce PD risk",
		subtitle = "ordered by FDR-adjusted p-value"
	) +
	plot_layout(
		widths = c(10, 3, 3),
		guides = 'collect'
	) & theme(legend.position = 'bottom')

ggsave(
	filename = "scanning_results_signif_decr_PD_risk.png", width = 10, height = 6)
