The directory contains scripts, data, and other relevant information on analysis
of the results of drugome-wise association study (DWAS) focusing on all the
ATC codes at level2; with Parkonson's disease (PD) as outcome.

The manuscript is on Overleaf.

## SCRIPTS

The scripts without numbers in the beginning were used inside the SAFE server,
on the raw data, to analyse and create descriptive tables. Scripts with the
numbers in the beginning were ran outside SAFE server, using only the results,
i.e., data available in the **DATA** folder.

### analysis.R

This is the main script. Cox regression with time-varying covariates, using
age of the individuals as time scale. Runs analyses on the entire dataset as
well as on the sex-stratified datasets.

### export_results_atc2level.R

The raw results saved by `analysis.R` were merged with help datasets:    
the ATC catalog and    
number of users per drug group.

This scripts produces the main results files, `expose_2_prescr_PD-4prescr-levo-mao-b_age-time-scale_results_atc2level_all.csv`, `expose_2_prescr_PD-4prescr-levo-mao-b_age-time-scale_results_atc2level_men.csv`, and `expose_2_prescr_PD-4prescr-levo-mao-b_age-time-scale_results_atc2level_women.csv`.

### check_descriptive_stats.R

Checking some information about the dataset, creating descriptive stats.

### checking_interaction.R

Run the interaction analyses for the chosen ATC groups.

### calc_n_users_per_drug_group.R

Calculate the exact number of users (both PD and non-PD) per ATC group.

### 01_coarse_analysis_pooled.R

Check results on analyses ran on the entire dataset. 
Plot results, create tables, and save for later.

### 02_coarse_analysis_sex-strat.R

Check results of analyses ran on the sex-stratified datasets.
Plot results, create tables, and save for later.


## RESULTS:

### `show_results_plotly.Rmd` and `show_results_plotly_sex_strata.Rmd`

Scripts that produce an interactive visualization of all the results.

Outputs:

- [show_results_plotly_atc_level2_pooled.html](show_results_plotly_atc_level2_pooled.html)
- [show_results_plotly_atc_level2_sex_strat.html](show_results_plotly_atc_level2_sex_strat.html)

## DATA:

Check the details in [data-folder README](DATA/00_README.md)
