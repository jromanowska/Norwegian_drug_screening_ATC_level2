The repository contains scripts, data, and other relevant information on analysis
of the results of drug-wise association study (DWAS) focusing on all the
ATC codes at level2; with Parkonson's disease (PD) as outcome.

The manuscript is published in Neurology:
"Association Between Use of Any of the Drugs Prescribed in Norway and the Subsequent Risk of Parkinson Disease:
A Drug-wide Association Study", Julia Romanowska, Kjetil Bjornevik, Marianna Cortese, Julia A Tuominen,
Magne Solheim, Asieh Abolpour Mofrad, Jannicke Igland, Clemens R Scherzer, Trond Riise, _Neurology_, **2023**,
[DOI: 10.1212/WNL.0000000000207899](https://n.neurology.org/content/early/2023/10/10/WNL.0000000000207899)

The raw files and results are available in the GitHub repo:
https://github.com/jromanowska/Norwegian_drug_screening_ATC_level2
and described below.

## SCRIPTS

The scripts in the `analysis_run` folder were used inside the SAFE server,
on the raw data, to analyse and create descriptive tables. Scripts with the
numbers in the beginning were ran outside SAFE server, using only the results,
i.e., data available in the **DATA** folder.

### `analysis_run` folder

- `data_preparation_db.R`    
Prepare dataset for Cox regression, all groups on ATC level2.

- `data_preparation_C09_level3_db.R`    
Prepare dataset for Cox regression, only the chosen sub-groups of C09.

- `data_preparation_dose-response.R`    
Prepare dataset for Cox regression, with number of prescriptions as a dose
proxy.

- `analysis.R`    
This is the script with main part of analysis procedure: Cox regression with
time-varying covariates, using age of the individuals as time scale. Other
scripts ending with _"analyses"_ set some parameters and call this one. 

- `main_analyses.R`    
Runs analyses on the entire dataset as well as on the sex-stratified
datasets. No time lag, using 2 prescriptions as exposure definition.

- `sensitivity_analyses.R`    
Runs analyses on the entire dataset as well as on the sex-stratified
datasets. No time lag, using 4 or 8 prescriptions as exposure definition.

- `time_lag_analyses.R`    
Runs analyses on the entire dataset as well as on the sex-stratified
datasets. Using 2 prescriptions as exposure definition, time lags of 5, 8,
and 10 years.

- `dose-response_analyses.R`    
Runs analyses for varying doses of drug group (using number of prescr.
as a proxy for dose).

- `level3_C09_analyses.R`    
Runs analyses on the two chosen sub-levels of C09 group.

- `export_results_atc2level.R`    
The raw results saved by `analysis.R` were merged with help datasets:    
the ATC catalog and    
number of users per drug group.

This scripts produces the main results files: 
- `expose_[n]_prescr_PD-4prescr-levo-mao-b_age-time-scale_results_atc2level_time-lag[x]years_all.csv`,
- `expose_[n]_prescr_PD-4prescr-levo-mao-b_age-time-scale_results_atc2level_time-lag[x]years_men.csv`, and
-`expose_[n]_prescr_PD-4prescr-levo-mao-b_age-time-scale_results_atc2level_time-lag[x]years_women.csv`.

Where _n_ is 2, 4, or 8; and _x_ is 0, 5, 8, or 10.

- `check_descriptive_stats.R`    
Checking some information about the dataset, creating descriptive stats.

- `checking_interaction.R`    
Run the interaction analyses for the chosen ATC groups.

- `calc_n_users_per_drug_group.R`    
Calculate the exact number of users (both PD and non-PD) per ATC group.

### `SCRIPTS` folder

- `01_coarse_analysis_pooled.R`    
Check results of the original analyses ran on the entire dataset. 
Plot results, create tables, and save for later.

- `02_coarse_analysis_sex-strat.R`    
Check results of the original analyses ran on the sex-stratified datasets.
Plot results, create tables, and save for later.

- `03_analysis_results_time-lag.R`    
Plot results, create tables, from time-lagged analyses.

- `04_analysis_results_sensitivity.R`    
Plot results, create tables, from sensitivity analyses.

- `05_dose-response_analyses.R`    
Plot results, create tables, from dose-response analyses.


## RESULTS:

### `show_results_plotly.Rmd` and `show_results_plotly_sex_strata.Rmd`

Scripts that produce an interactive visualization of all the results.

Outputs:

- [show_results_plotly_atc_level2_pooled.html](RESULTS/show_results_plotly_atc2level.html)
- [show_results_plotly_atc_level2_sex_strat.html](RESULTS/show_results_plotly_ATC2level_sex_strata.html)

### `nice_results_atc_group3_C09_*.rds`

R-data with all the results from the sub-group analysis.

### Dose-response analysis

All results can be [viewed here](RESULTS/dose_response_table_with_plots.html).

## DATA:

Check the details in [data-folder README](DATA/)
