# DESCRIPTION: Creation of dataset; based on Magne's script
#    S:/Project/DRONE/Magne/test4/Data_prep_4.Rmd
#	   Ta ut alle med minst en N04-resept i 2004
#    Minst 10 PD-caser per medikament som skal være med i screening
#    **PD-definisjon:**
#	    * minst 4 Levodopa-resepter + minst en resept (av hvilket som helst legemiddel)
#       med refusjonskode G20 /pktNr 16, eller
#     * minst 4 MAO-B-resepter
#
#    **PD-dato:** dato for første N04 resept
#
#   HERE: using the database implementation
#   IMPORTANT UPDATE: dynamically changing the min.number of prescriptions
#                     required for exposure
#
# AUTHOR: Julia Romanowska
# DATE CREATED: 2023-01-31
# DATE LAST MODIFIED: 2023-06-01

# SETUP ----
library(tidyverse)
library(lubridate)
library(data.table)
library(here)
library(RSQLite)
library(dbplyr)

all_prescriptions_db <- file.path("..", "..", "DrugScreeningDB.sqlite")
prescriptions_by_ATC_level2_file <- here("DATA",
																				 "prescriptions_by_ATC_level2_first10.rds")
prescriptions_exposure <- 8

ssb_censor_date_file <- here("DATA", "SSB_censor_date_event.rds")
ssb_edu_level_file <- here("DATA", "SSB_edu_level.rds")
PD_patients_identified_file <- here("DATA", "PD_patients_identified.rds")

# READ DATA -----------------
db_connection <- DBI::dbConnect(RSQLite::SQLite(), all_prescriptions_db)

dbListTables(db_connection)

## SSB
ssb_npr_data <- readRDS(
	file.path("..", "..", "Copy_original", "DATA_CLEANING",
						"SSB-NPR_data_cleaned.rds")
)
names(ssb_npr_data)

ssb_birth_death <- ssb_npr_data$PDB2687_W19_1693_BEFOLKNING_FASTE_202005.dta
ssb_birth_death
ssb_birth_death <- ssb_birth_death %>%
	select(-fraland_aar, -tilland_aar, -innalder)

setDT(ssb_birth_death)

ssb_regstat <- ssb_npr_data$PDB2687_W19_1693_BEFOLKNING_REGSTAT_202005.dta
table(ssb_regstat$regstat_04)
ssb_registered2004 <- ssb_regstat$pasientlopenr_pdb2687

# registration status
if(!file.exists(ssb_censor_date_file)){
	# need to check when each person emmigrated - we will censor them at this date
	#  regstat == 1 - living in Norway
	#  regstat == 3 - emmigrated
	ssb_censor_date <- ssb_regstat %>%
		select(PasientLopeNr = pasientlopenr_pdb2687, everything()) %>%
		rowwise(PasientLopeNr) %>%
		# paste all the statuses together into one string
		summarise(
			regstat = paste0(
				as.character(c_across(regstat_04:regstat_18)),
				collapse = "")
		) %>%
		# find first change of status and set it to censored date
		mutate(
			censor_year = ifelse(
				!is.na(str_locate(regstat, "3")[,1]),
				yes = 2003 + str_locate(regstat, "3")[,1] - 1,
				no = ifelse(
					!is.na(str_locate(regstat, "5")[,1]),
					yes = 2003 + str_locate(regstat, "5")[,1] - 1,
					no = 2019
				)
			),
			censor_event = ifelse(
				!is.na(str_locate(regstat, "3")[,1]),
				yes = "emmigrated",
				no = ifelse(
					!is.na(str_locate(regstat, "5")[,1]),
					yes = "died",
					no = "alive"
				)
			)
		) %>% ungroup()
	ssb_censor_date
	ssb_censor_date %>%
		count(censor_event)
	
	saveRDS(ssb_censor_date, file = ssb_censor_date_file)
} else {
	ssb_censor_date <- readRDS(ssb_censor_date_file)
}

rm(ssb_regstat)

# education
if(!file.exists(ssb_edu_level_file)){
	ssb_edu <- ssb_npr_data$PDB2687_W19_1693_UTDANNING_BU_202005.dta
	ssb_edu
	ssb_edu <- ssb_edu %>%
		select(bu_nivaa_2004, bu_gruppe_2004, PasientLopeNr = pasientlopenr_pdb2687)
	ssb_edu
	
	# we don't need so many education levels
	ssb_edu_key <- ssb_edu %>%
		distinct(bu_nivaa_2004, bu_gruppe_2004)
	ssb_edu_key
	ssb_edu %>%
		count(bu_gruppe_2004)
	ssb_edu <- ssb_edu %>%
		mutate(
			edu_level = case_when(
				(is.na(bu_nivaa_2004) | bu_nivaa_2004 == 0) ~ 0,
				bu_nivaa_2004 %in% c(1, 2) ~ 1,
				bu_nivaa_2004 %in% c(3, 4, 5) ~ 2,
				bu_nivaa_2004 == 6 ~ 3,
				bu_nivaa_2004 %in% c(7, 8) ~ 4
			)
		)
	ssb_edu
	ssb_edu %>%
		count(edu_level)
	
	saveRDS(ssb_edu, file = ssb_edu_level_file)
	
} else {
	ssb_edu <- readRDS(ssb_edu_level_file)
}
rm(ssb_npr_data)

# IDENTIFYING PD PATIENTS ----
prescriptions <- tbl(db_connection, "Prescription_data")

if(!file.exists(PD_patients_identified_file)){
	# we will gather all N04BA and all N04BD into two categories
	#   because we don't differentiate those
	prescriptions_N04 <- prescriptions %>%
		filter(ATC_code %like% "N04%") %>%
		collect()
	
	setDT(prescriptions_N04)
	
	prescriptions_N04[ATC_code %in% c("N04BA02", "N04BA03"), ATC_code := "N04BA"]
	prescriptions_N04[ATC_code %in% c("N04BD01", "N04BD02", "N04BD03"), ATC_code := "N04BD"]
	
	setkey(prescriptions_N04, ATC_code)
	prescriptions_N04[, UtleveringsDato := ymd(
		paste(prescr_year, prescr_mnth, "15", sep = "-")
	)]
	prescriptions_N04
	
	# calculate cummulative number of prescriptions
	setorder(prescriptions_N04, indiv_ID, ATC_code, UtleveringsDato)
	prescriptions_N04[, kum_resept := seq(.N), 
							 keyby = .(indiv_ID, ATC_code)]
	
	# CRITERIUM 1: min.4 prescr.of levo and PD reimbursement code on any prescription
	patients_crit1a <- prescriptions %>%
		filter(refusjonICD == "G20" | refusjonPktNr == 16) %>%
		distinct(indiv_ID) %>%
		collect() %>%
		pull(indiv_ID)
	
	patients_crit1b <- prescriptions_N04[
		(kum_resept >= 4 & ATC_code == "N04BA"),
		.SD[1],
		by = "indiv_ID"
	]$indiv_ID
	patients_crit1 <- intersect(
		patients_crit1a,
		patients_crit1b
	)
	
	length(patients_crit1)
	
	# CRITERIUM 2: min.4 prescr.of any MAO-B inhibitors
	patients_crit2 <- prescriptions_N04[
		(kum_resept >= 4 & ATC_code == "N04BD"),
		.SD[1],
		by = "indiv_ID"
	]$indiv_ID
	length(patients_crit2)
	
	# COMBINE - any of the 2 criteria
	patients_crit_any <- union(patients_crit1, patients_crit2)
	length(patients_crit_any)
	
	# patients_PD <- data.table(
	# 	indiv_ID = patients_crit_any,
	# 	park_yn = TRUE
	# )
	
	# find first date for any N04 prescription
	patients_PD_first_date <- prescriptions_N04[
		(str_sub(ATC_code, 1, 3) == "N04") & (indiv_ID %in% patients_crit_any),
		.(first_N04_date = min(UtleveringsDato)),
		by = .(indiv_ID)
	]
	
	# CHECKING
	# how many patients that had their first N04 prescr. after 2004?
	uniqueN(patients_PD_first_date[
		first_N04_date >= ymd('2005-01-01')
	])
	
	saveRDS(patients_PD_first_date, file = PD_patients_identified_file)
} else {
	patients_PD_first_date <- readRDS(PD_patients_identified_file)
}

## GROUP BY ATC LEVEL 2 ----
if(!file.exists(prescriptions_by_ATC_level2_file)){
	prescriptions_with_level2 <- prescriptions %>%
		mutate(ATC_level2 = substr(ATC_code, 1, 3)) %>%
		count(indiv_ID, ATC_level2) %>%
		collect()

# we consider exposure only when there are min.2, 4, or 8 prescriptions
# then, we need the date of first as the exposure date
	setDT(prescriptions_with_level2)
	
	individuals_exposed <- prescriptions_with_level2[
		n >= prescriptions_exposure,
	]
	prescriptions_with_level2_exposed <- prescriptions %>%
		select(indiv_ID, prescr_year, prescr_mnth, ATC_code) %>%
		mutate(ATC_level2 = substr(ATC_code, 1, 3)) %>%
		group_by(indiv_ID, ATC_level2) %>%
		arrange(indiv_ID, ATC_level2, prescr_year, prescr_mnth) %>%
		mutate(prescr_date = paste(prescr_year, prescr_mnth, "15", sep = "-")) %>%
		group_by(prescr_date) %>%
		mutate(prescr_no = row_number()) %>%
		filter(prescr_no == 1) %>% # this will take only first prescription each month
		ungroup() %>%
		mutate(prescr_no = row_number()) %>%
		filter(prescr_no == prescriptions_exposure) %>% # this will take only the rows with the actual prescription that gives exposed/non-exposed category
		collect()
	
	prescriptions_with_level2_exposed
	setDT(prescriptions_with_level2_exposed)

	prescriptions_with_level2_exposed[
		, start_exposure_date := ymd(
			paste(prescr_year, prescr_mnth, "15", sep = "-")
		)
	]
} else {
	# this file was saved in 'save_prescriptions_ATC_level2_first10_prescr.R'
	prescriptions_with_level2 <- readRDS(prescriptions_by_ATC_level2_file)
	
	# I need to merge first the prescriptions that are from the same month!
	prescriptions_with_level2[
		, kum_resept_per_mnth := seq(.N), 
		keyby = .(PasientLopeNr, ATCKode2, UtleveringsDato)
	]
	# filter by taking only first prescription of the group for each month
	first_prescriptions_with_level2 <- prescriptions_with_level2[
		kum_resept_per_mnth == 1,
	][
		, kum_resept := seq(.N),
		keyby = .(PasientLopeNr, ATCKode2)
	]
	first_prescriptions_with_level2
	# just checking:
	first_prescriptions_with_level2[, .N, by = kum_resept_per_mnth][, .N, by = N]
	
	prescriptions_with_level2_exposed <- first_prescriptions_with_level2[
		kum_resept == prescriptions_exposure,
		.(indiv_ID = PasientLopeNr, n = kum_resept, ATC_code = ATCKode,
			ATC_level2 = ATCKode2, start_exposure_date = UtleveringsDato)
	]
}

## CREATE EXPOSURE VARIABLE -----------
# time to exposure = time from the start of observation
start_date <- lubridate::ymd("2005-01-01")
prescriptions_with_level2_exposed[,
	med_time := as.numeric(start_exposure_date - ..start_date)
]

# join with some other info about individuals
dataset_long <- merge.data.table(
	prescriptions_with_level2_exposed,
	patients_PD_first_date,
	by = "indiv_ID",
	all = TRUE
)
dataset_long[, park_yn := !is.na(first_N04_date)]
dataset_long

# count (roughly) number of PD-patients using each drug
PD_patients_per_drug_group <- dataset_long[
	park_yn == TRUE,
	.N,
	by = ATC_level2
]
PD_patients_per_drug_group[
	N >= 5
][
	order(N)
]

retain_ATC_level2 <- PD_patients_per_drug_group[
	N >= 5, ATC_level2
]
retain_ATC_level2 <- str_subset(retain_ATC_level2, "N04", negate = TRUE)

## RESTRUCTURE DATASET: now, only one row per individual ------
data_restr <- dcast(
	dataset_long,
	indiv_ID + park_yn + first_N04_date ~ ATC_level2,
	value.var = c("med_time")
)
dim(data_restr)
names(data_restr)

data_restr[,1:10]
names_not_atc <- names(data_restr)[1:3]
names_retain <- c(names_not_atc, sort(retain_ATC_level2))
names_retain

data_restr <- data_restr[
	, ..names_retain
]
data_restr[,1:10]

# check - PD status
data_restr[, .N, by = park_yn]

## REMOVE PD PATIENTS WHO HAD N04 PRESCR.IN 2014
# NOTE: cannot remove here because it will cause errors when adding these
#   persons afterwards from the SSB dataset!
## CHECK WHICH PD PATIENTS HAD N04 PRESCR.IN 2014 ----
data_restr[!is.na(first_N04_date), 1:10]

data_restr <- data_restr[
	!is.na(first_N04_date),
	first_N04_before_start := if_else(
		first_N04_date < ymd("2005-01-01"),
		TRUE,
		FALSE
	)
]
data_restr <- data_restr[
	is.na(first_N04_date),
	first_N04_before_start := FALSE
]
dim(data_restr)
data_restr[!is.na(first_N04_date), 1:10]

data_restr[, .N, by = .(park_yn, first_N04_before_start)]

# re-calculate PD onset to days after start of follow-up
data_restr[,
	 time_to_onset := ifelse(
		 	park_yn & !first_N04_before_start,
	 		yes = as.numeric(first_N04_date - ..start_date),
		 	no = NA
	 )
]
data_restr[!is.na(first_N04_date)]

## ADDING DATA FROM SSB ----
# merging files - there are some people that are not in reseptregisteret, but
#  have data in SSB
data_restr_ssb <- merge(
	data_restr,
	ssb_birth_death[
		, .(foedselsaar, doeds_aar, sex = kjoenn, indiv_ID = pasientlopenr_pdb2687)],
	by = "indiv_ID",
	all = TRUE
)
dim(data_restr_ssb)
names(data_restr_ssb)

data_restr_ssb[, .N, by = park_yn]

# all these with park_yn == NA are not PD patients
data_restr_ssb[, park_yn := if_else(
	is.na(park_yn),
	FALSE,
	park_yn
)]

# there are some persons that are not in SSB - these were probably not registered
#  as living in Norway on 1.1.2004
sum(is.na(data_restr_ssb$foedselsaar))
data_restr_ssb[is.na(foedselsaar), .N, by = .(park_yn, first_N04_before_start)]

data_restr_ssb[is.na(foedselsaar) & indiv_ID %in% ssb_registered2004]

data_restr_ssb <- data_restr_ssb[!is.na(foedselsaar)]

## REMOVE PEOPLE WHO HAD ANY N04 in 2004! ----
data_restr_ssb[, .N, by = first_N04_before_start]
data_restr_ssb <- data_restr_ssb[
	is.na(first_N04_before_start) | !first_N04_before_start,
]
data_restr_ssb[, .N, by = first_N04_before_start]

# get the death date from Prescription registry for those who are there
death_dates_NorPD <- tbl(db_connection, "Individuals_data") %>%
	select(indiv_ID, death_year, death_mnth) %>%
	collect()
death_dates_NorPD %>% count(death_mnth)
death_dates_NorPD %>% filter(is.na(death_year)) %>% count(death_mnth)
# I think I've accidentally added '07' to death month
# for _everyone_, not only those who died!
# need to recode that '07'!
setDT(death_dates_NorPD)

data_restr_ssb <- merge(
	data_restr_ssb,
	death_dates_NorPD,
	by = "indiv_ID",
	all.x = TRUE
)
data_restr_ssb[is.na(death_year), .N, by = death_mnth]
data_restr_ssb[is.na(death_year) & death_mnth == "07", .N, by = park_yn]
data_restr_ssb[is.na(death_year) & death_mnth != "07", .N, by = park_yn]

data_restr_ssb[, .N, by = first_N04_before_start]

# age in 2004
data_restr_ssb[, age := 2004 - foedselsaar]
summary(data_restr_ssb$age)
nrow(data_restr_ssb[age > 25])

# adding education
ssb_edu_dt <- as.data.table(
	ssb_edu %>%
		select(indiv_ID = PasientLopeNr, edu_level)
)

data_restr_ssb <- merge(
	data_restr_ssb,
	ssb_edu_dt,
	by = "indiv_ID",
	all.x = TRUE
)
data_restr_ssb[, .N, by = edu_level][order(N)]
data_restr_ssb[, .N, by = .(park_yn, edu_level)][order(park_yn, N)]
data_restr_ssb[, .N, by = .(park_yn, first_N04_before_start)]

## SETTING STATUS TIME -----
# censoring:
#   - those who are alive will get 2020.1.1
#   - those who died, will get their death date
#   - those who emmigrated will get the 1.7 of the emmigration year
data_restr_ssb <- merge(
	data_restr_ssb,
	setDT(ssb_censor_date %>%
					select(-regstat) %>%
					rename(indiv_ID = PasientLopeNr)),
	by = "indiv_ID",
	all.x = TRUE
)
data_restr_ssb[,
	cens_date := ymd("2020-01-01")
]
data_restr_ssb[
	censor_event == "died",
	cens_date := ymd(paste(censor_year, death_mnth, "15", sep = "-"))
]
# NB! important with order here: if emmigration is before death,
#     then this should be censoring date!
data_restr_ssb[
	censor_event == "emmigrated",
	cens_date := ymd(paste(censor_year, "07", "15", sep = "-"))
]

data_restr_ssb[, .N, by = park_yn]
data_restr_ssb[, .N, by = .(park_yn, first_N04_before_start)]
data_restr_ssb[is.na(censor_year)]
data_restr_ssb[is.na(death_year), .N, by = censor_event]

data_restr_ssb <- data_restr_ssb[, -c("doeds_aar", "death_year", "death_mnth", "first_N04_before_start")]

## TIME UNDER RISK -----
# for those who had PD: time = PD onset date - start
# for others: time = status date - start

data_restr_ssb[
	, status_date := pmin(cens_date, first_N04_date, na.rm = T)
]

data_restr_ssb[, time_risk := status_date - ..start_date]
data_restr_ssb

## READY FOR ANALYSIS - WRITE THE DATASET TO THE DISK --------
saveRDS(
	data_restr_ssb,
	file = here(
		"DATA",
		paste0("dataset_ready_for_analysis_exposure", prescriptions_exposure,".rds")
	)
)


# DISCONNECT ----
dbDisconnect(db_connection)

