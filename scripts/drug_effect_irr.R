#------------------------------------------------------------------------------------------------------
# INTER-RATER AGREEMENT: ICC AT 90 MIN
#------------------------------------------------------------------------------------------------------

library(irr)

icc_data <- monrat %>%
  filter(!(study_id == '1405' & sub == 'NYU')) %>%
  filter(timepoint == 90) %>%
  mutate(
    SessMonRat_OverallDrugEffect = as.numeric(SessMonRat_OverallDrugEffect),
    SessMonRat_RaterCode         = as.numeric(SessMonRat_RaterCode)
  ) %>%
  filter(SessMonRat_RaterCode %in% c(1, 2)) %>%
  select(study_id, ID, session, SessMonRat_RaterCode, SessMonRat_OverallDrugEffect) %>%
  pivot_wider(
    names_from  = SessMonRat_RaterCode,
    values_from = SessMonRat_OverallDrugEffect,
    names_prefix = "rater_"
  ) %>%
  filter(!is.na(rater_1), !is.na(rater_2))

# ICC: two-way mixed, absolute agreement, single measures
icc_result <- icc(
  icc_data %>% select(rater_1, rater_2),
  model = "twoway",
  type  = "agreement",
  unit  = "single"
)

print(icc_result)

# Also compute raw agreement and kappa
cat("\nExact agreement: ", 
    round(mean(icc_data$rater_1 == icc_data$rater_2), 3), "\n")
cat("Within-1 agreement: ",
    round(mean(abs(icc_data$rater_1 - icc_data$rater_2) <= 1), 3), "\n")

# Repeat across all timepoints for completeness
icc_by_tp <- map_dfr(c(30, 60, 90, 120), function(tp) {
  d <- monrat %>%
    filter(!(study_id == '1405' & sub == 'NYU')) %>%
    filter(timepoint == tp) %>%
    mutate(
      SessMonRat_OverallDrugEffect = as.numeric(SessMonRat_OverallDrugEffect),
      SessMonRat_RaterCode         = as.numeric(SessMonRat_RaterCode)
    ) %>%
    filter(SessMonRat_RaterCode %in% c(1, 2)) %>%
    select(study_id, ID, session, SessMonRat_RaterCode, SessMonRat_OverallDrugEffect) %>%
    pivot_wider(
      names_from   = SessMonRat_RaterCode,
      values_from  = SessMonRat_OverallDrugEffect,
      names_prefix = "rater_"
    ) %>%
    filter(!is.na(rater_1), !is.na(rater_2))
  
  res <- icc(d %>% select(rater_1, rater_2),
             model = "twoway", type = "agreement", unit = "single")
  
  tibble(
    timepoint      = paste0(tp, " min"),
    n_pairs        = nrow(d),
    icc            = round(res$value, 3),
    icc_lo         = round(res$lbound, 3),
    icc_hi         = round(res$ubound, 3),
    exact_agreement = round(mean(d$rater_1 == d$rater_2), 3),
    within_1        = round(mean(abs(d$rater_1 - d$rater_2) <= 1), 3)
  )
})

print(icc_by_tp)
