## Session Monitor Rating Data Combined
## nsepeda1@jh.edu

###*****************************************************************************************************
### INSTALL AND LOAD REQUIRED PACKAGES
###*****************************************************************************************************
required.packages <- c("dplyr","qualtRics","data.table","tidyverse","officer","openxlsx","httr") #list packages used in script
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])] #checks for new packages from required.packages
if(length(new.packages)) install.packages(new.packages) # installs any packages don't yet have
lapply(required.packages, require, character.only = TRUE) #load specified packages


###*****************************************************************************************************
### SETUP AND READ DATA
###*****************************************************************************************************

## Clear environment
start_fresh <- 1
if (start_fresh == 1){
  rm(list=ls())
}

## Set working directory
setwd('/Users/nsepeda/Documents/00xx_DHProjects/DrugEffectPrediction')

## Read combined session monitor rating data (studies 0014 - 1909)
monrat <- read_csv('Combined_MonitorSessionRating.csv')
guidefx <- read_csv('GuideEffectsData.csv')


###*****************************************************************************************************
### CLEAN DATA SET
###*****************************************************************************************************
#------------------------------------------------------------------------------------------------------
# Guide Effects Data
#------------------------------------------------------------------------------------------------------
guidefx_clean <- guidefx %>%
  filter(guide_type == 1) %>%
  rename(session = timepoint)

#------------------------------------------------------------------------------------------------------
# Within-Session Monitor Rating Data (Combined from studies 0014 to 1909)
#------------------------------------------------------------------------------------------------------
monrat_clean <- monrat %>%
  filter(!(study_id == '1405' & sub == 'NYU')) %>% #remove NYU 1405 vols
  select(-c(sub, substance, condition, dose, code)) %>%
  mutate(study_name = ifelse(study_name == "Beginner Mediatator Psilocybin", 
                             "Beginner Meditator Psilocybin", 
                             study_name)) %>%
  mutate(study_name = ifelse(study_name == "Religious Professionals and Psiolocybin", 
                             "Religious Professionals and Psilocybin", 
                             study_name)) %>%
  mutate(session = ifelse(study_id == 1610, 1, session))

write.csv(monrat_clean, 'Combined_MonitorSessionRating_CLEAN.csv', row.names = F)


## Calculate Maximum drug effect rating between raters
monrat_max <- monrat_clean %>%
  select("study_id", "study_name", "ID", "session", "date", "time", "timepoint", 
         "SessMonRat_RaterCode", "SessMonRat_initials", "SessMonRat_OverallDrugEffect") %>%
  mutate(SessMonRat_OverallDrugEffect = as.numeric(SessMonRat_OverallDrugEffect)) %>%
  group_by(study_id, study_name, ID, session, timepoint) %>%
  summarise(Max_OverallDrugEffect = suppressWarnings(max(SessMonRat_OverallDrugEffect, na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(Max_OverallDrugEffect = ifelse(is.infinite(Max_OverallDrugEffect), NA_real_,
                                        Max_OverallDrugEffect)) %>%
  mutate(timepoint = as.numeric(timepoint),
         study_id = as.numeric(study_id),
         ID = as.numeric(ID),
         session = as.numeric(session),
         Max_OverallDrugEffect = as.numeric(Max_OverallDrugEffect)) %>%
  arrange(study_id, ID, session, timepoint)

write.csv(monrat_max, 'Combined_MonitorSessionRating_MAX.csv', row.names = F)


## Calculate Average drug effect rating between raters
monrat_avg <- monrat_clean %>%
  select("study_id", "study_name", "ID", "session", "date", "time", "timepoint",
         "SessMonRat_RaterCode", "SessMonRat_initials", "SessMonRat_OverallDrugEffect") %>%
  mutate(SessMonRat_OverallDrugEffect = as.numeric(SessMonRat_OverallDrugEffect)) %>%
  group_by(study_id, study_name, ID, session, timepoint) %>%
  summarise(Avg_OverallDrugEffect = mean(SessMonRat_OverallDrugEffect, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(Avg_OverallDrugEffect = ifelse(is.nan(Avg_OverallDrugEffect), NA_real_,
                                        Avg_OverallDrugEffect)) %>%
  mutate(timepoint = as.numeric(timepoint),
         study_id = as.numeric(study_id),
         ID = as.numeric(ID),
         session = as.numeric(session),
         Avg_OverallDrugEffect = as.numeric(Avg_OverallDrugEffect)) %>%
  arrange(study_id, ID, session, timepoint)
  
write.csv(monrat_avg, 'Combined_MonitorSessionRating_AVG.csv', row.names = F)


#------------------------------------------------------------------------------------------------------
# Assess Similar vs. Unique Combinations of Vol Data Between Guide Effects and Mon Rat Data
#------------------------------------------------------------------------------------------------------
## pull just a single timepoint from monrat_max to assess overlap with guidefx data
monrat_max_vol <- monrat_max %>% filter(timepoint == 30)

## Assess overlap vs uniqueness of study_id, study_name, ID, and session combinations

# Define grouping keys
keys <- c("study_id", "study_name", "ID", "session")

# Get distinct key combos from each dataset
guidefx_keys <- guidefx_clean %>% distinct(across(all_of(keys)))
monrat_keys <- monrat_max_vol %>% distinct(across(all_of(keys)))

# Overlap: in both
overlap <- inner_join(guidefx_keys, monrat_keys, by = keys)

# Unique to guidefx_clean
only_guidefx <- anti_join(guidefx_keys, monrat_keys, by = keys)

# Unique to monrat_max_vol
only_monrat <- anti_join(monrat_keys, guidefx_keys, by = keys)

# Summary
cat("Overlapping combinations:        ", nrow(overlap), "\n")
cat("Unique to guidefx_clean:         ", nrow(only_guidefx), "\n")
cat("Unique to monrat_max_vol:        ", nrow(only_monrat), "\n")


#------------------------------------------------------------------------------------------------------
# Assess Unique Vol ID in Comparing monrat and guidefx data sets
#------------------------------------------------------------------------------------------------------
id_comparison <- data.frame()

for (i in unique(monrat_max$study_id)){
  ids_monrat <- monrat_max %>% filter(study_id == i) %>% pull(ID) %>% unique()
  ids_guidefx <- guidefx_clean %>% filter(study_id == i) %>% pull(ID) %>% unique()
  
  only_monrat <- setdiff(ids_monrat, ids_guidefx)
  only_guidefx <- setdiff(ids_guidefx, ids_monrat)
  
  id_comparison <- bind_rows(id_comparison, data.frame(
    study_id = i,
    unique_id_monrat = ifelse(length(only_monrat) == 0, NA_character_, paste(only_monrat, collapse = ", ")),
    unique_id_guidefx = ifelse(length(only_guidefx) == 0, NA_character_, paste(only_guidefx, collapse = ", "))
  ))
}

