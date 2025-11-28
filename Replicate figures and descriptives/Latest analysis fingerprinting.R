# ------------------------------------------------------------------------
# Author: Martina Occelli (mo386@cornell.edu)
# Project: Whose bias? DNA fingerprinting paper
# Description: 03 - Descriptive manuscript (October 2025)
# ------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(plyr)
library(haven) # read dta files
library(readr) # read excel/csv files
library(readxl)
library(writexl)
library(dplyr)
library(ggplot2) # visualization
library(miceadds) #clustered errors
library(estimatr) #clustered errors
library(DescTools) #winsorizing
library(scales)
library(gtsummary) # to use tbl_summary for compare populations
library(lmtest)
library(sandwich)
library(dotwhisker)

# Load datasets
main_survey_clean_Sep9_fingerprinting <- read_csv("~/Desktop/main_survey_clean_Sep9-fingerprinting.csv")
main_data <- main_survey_clean_Sep9_fingerprinting

main_data <- main_data %>%
  select(2:91)

fingerprint_data <- read_excel("~/Google Drive/My Drive/!!CORNELL/ILCI /Priority Setting/DNA fingerprinting_Project/!Manuscripts/
                               !Whose bias - econ finger manuscript /Sergio Data-Analysis/fingerprint_updated_data.xlsx")

# Identify unique IDs in both datasets
unique_main_ids <- unique(main_data$farm_id)
unique_fingerprint_codes <- unique(fingerprint_data$farm_code)

ids_not_in_fingerprint <- setdiff(unique_main_ids, unique_fingerprint_codes) # there are 26 HH for which we have no fingerprinting results

# Filter main survey data excluding ids which are not in the fingerprinting main dataset
main_data_matched <- main_data[(main_data$farm_id %in% fingerprint_data$farm_code), ] 
unique_matched_ids <- unique(main_data_matched$farm_id) # double check that we have 668 unique observations

# Table A3, sample summary statistics (survey only)
table(main_data_matched$treatment_list)
table(fingerprint_data$treated)

table <- aggregate(farm_code ~ treated, data = fingerprint_data, FUN = function(x) length(unique(x))) # how many treated and control we have in fingerprinting
print(table)

table_2 <- aggregate(farm_id ~ treatment_list, data = main_data_matched, FUN = function(x) length(unique(x))) # how many treated and control we have in survey
print(table_2)

# The two don't match (we have 7HH shifting from treatment in fingerprinting to control in survey)
# Why? Who are they?

fp_treatment_map <- unique(fingerprint_data[, c("farm_code", "treated")]) #fingerprinting data
main_treatment_map <- unique(main_data_matched[, c("farm_id", "treatment_list")]) # survey data

names(fp_treatment_map) <- c("farm_id", "treated_fp")
names(main_treatment_map) <- c("farm_id", "treated_main")

merged_treatment <- merge(fp_treatment_map, main_treatment_map, by = "farm_id", all = TRUE)

# Find rows where treatment assignments don't match
mismatched_ids <- merged_treatment[merged_treatment$treated_fp != merged_treatment$treated_main, ]

# Or where treatment is missing from one dataset
missing_treatment <- merged_treatment[is.na(merged_treatment$treated_fp) | is.na(merged_treatment$treated_main), ]

# Table by treatment based on the survey data
str(main_data_matched$tricot_check)
main_data_matched$tricot_check_recoded <- ifelse(main_data_matched$tricot_check == "1", 1,
                                                  ifelse(main_data_matched$tricot_check == "0", 0, NA))

t.test(main_data_matched$tricot_check_recoded ~ treatment_list, data = main_data_matched)

aggregate(main_data_matched$tricot_check_recoded ~ treatment_list, data = main_data_matched, 
          FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                              sd = sd(x, na.rm = TRUE)))


fingerprint_data$gender_recoded <- ifelse(fingerprint_data$gender == 2, 1,
                                          ifelse(fingerprint_data$gender == 1, 0, NA))
t.test(education ~ treated, data = fingerprint_data)

aggregate(education ~ treated, data = fingerprint_data, 
          FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                              sd = sd(x, na.rm = TRUE)))

# export main datasets
write.csv(fingerprint_data, 'fingerprint_data.csv')
write.csv(main_data_matched, 'survey_data_matched.csv')
write.csv(mismatched_ids, 'mismatch_fing_survey_treatment.csv')


