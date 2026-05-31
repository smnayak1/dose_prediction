## Drug Effect → Redosing Decision Rule
## Monitor-rated drug effect at 90 min predicting missed CME (MEQ < threshold)
## Decision rule: redose if avg drug effect ≤ 1 at 90 min
# When to Redose? An Empirical Decision Rule for In-Session Psilocybin Redosing  ;Based on Real-Time Session Monitoring



library(here)
library(tidyverse)
library(pROC)

#------------------------------------------------------------------------------------------------------
# READ DATA
#------------------------------------------------------------------------------------------------------

monrat_avg <- read_csv(here("data", "Combined_MonitorSessionRating_AVG.csv"))
guidefx    <- read_csv(here("data", "GuideEffectsData_CLEAN.csv"))

#------------------------------------------------------------------------------------------------------
# PREP
#------------------------------------------------------------------------------------------------------

cme_cutoff <- 0.60  # MEQ % max threshold — update as needed

guidefx_clean <- guidefx %>%
  select(study_id, study_name, ID, session, meq.total.pctmax) %>%
  mutate(across(c(study_id, ID, session), as.numeric))

df <- monrat_avg %>%
  filter(timepoint %in% c(30, 60, 90, 120)) %>%
  select(study_id, study_name, ID, session, timepoint, Avg_OverallDrugEffect) %>%
  inner_join(guidefx_clean, by = c("study_id", "study_name", "ID", "session")) %>%
  filter(!is.na(meq.total.pctmax), !is.na(Avg_OverallDrugEffect)) %>%
  mutate(across(c(study_id, ID, session, timepoint), as.numeric))

cat("Sessions:     ", df %>% distinct(study_id, ID, session) %>% nrow(), "\n")
cat("Participants: ", df %>% distinct(ID) %>% nrow(), "\n")
df %>% count(timepoint)

#------------------------------------------------------------------------------------------------------
# EMPIRICAL P(CME) BY DRUG EFFECT RATING AND TIMEPOINT
#------------------------------------------------------------------------------------------------------

tp_colors <- c("30" = "#2C6E8A", "60" = "#E05A2B",
               "90" = "#2CA87F", "120" = "#9B59B6")

emp_prob <- df %>%
  mutate(
    de_bin = round(Avg_OverallDrugEffect * 2) / 2,
    cme    = as.integer(meq.total.pctmax >= cme_cutoff)
  ) %>%
  group_by(timepoint, de_bin) %>%
  summarise(p_cme = mean(cme, na.rm = TRUE), n = n(), .groups = "drop") %>%
  filter(n >= 3)

ggplot(emp_prob, aes(x = de_bin, y = p_cme,
                     color = factor(timepoint),
                     group = factor(timepoint))) +
  geom_line(linewidth = 1.2) +
  geom_point(aes(size = n)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "#AAAAAA") +
  scale_color_manual(values = tp_colors,
                     labels = paste0(c(30, 60, 90, 120), " min"),
                     name   = "Timepoint") +
  scale_size_continuous(range = c(2, 7), name = "N sessions") +
  scale_x_continuous(breaks = seq(0, 4, by = 0.5)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(
    title = paste0("Empirical P(CME) by Drug Effect Rating  |  MEQ ≥ ", round(cme_cutoff * 100), "%"),
    x     = "Avg drug effect rating",
    y     = paste0("P(MEQ ≥ ", round(cme_cutoff * 100), "%)")
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank())

#------------------------------------------------------------------------------------------------------
# DECISION RULE SWEEP
# Positive = missed CME; flagged = drug effect <= threshold → redose
#------------------------------------------------------------------------------------------------------

redose_table <- map_dfr(c(30, 60, 90, 120), function(tp) {
  map_dfr(seq(0, 3, by = 0.5), function(thresh) {
    
    d <- df %>%
      filter(timepoint == tp) %>%
      group_by(study_id, ID, session, meq.total.pctmax) %>%
      summarise(drug_effect = mean(Avg_OverallDrugEffect, na.rm = TRUE), .groups = "drop") %>%
      filter(!is.na(drug_effect)) %>%
      mutate(
        missed_cme = as.integer(meq.total.pctmax < cme_cutoff),
        flagged    = as.integer(drug_effect <= thresh)
      )
    
    tp_count <- sum(d$missed_cme == 1 & d$flagged == 1)
    tn_count <- sum(d$missed_cme == 0 & d$flagged == 0)
    fp_count <- sum(d$missed_cme == 0 & d$flagged == 1)
    fn_count <- sum(d$missed_cme == 1 & d$flagged == 0)
    n_missed <- sum(d$missed_cme)
    n_hit    <- nrow(d) - n_missed
    
    tibble(
      timepoint    = paste0(tp, " min"),
      de_thresh    = thresh,
      n_sessions   = nrow(d),
      n_missed_cme = n_missed,
      n_flagged    = tp_count + fp_count,
      sensitivity  = round(tp_count / n_missed,              3),
      specificity  = round(tn_count / n_hit,                 3),
      ppv          = round(tp_count / (tp_count + fp_count), 3),
      npv          = round(tn_count / (tn_count + fn_count), 3),
      fpr          = round(fp_count / n_hit,                 3),
      fnr          = round(fn_count / n_missed,              3)
    )
  })
}) %>%
  mutate(across(c(sensitivity, specificity, ppv, npv, fpr, fnr),
                ~ ifelse(is.nan(.), NA_real_, .)))

print(redose_table, n = Inf)

# Heatmap
redose_table %>%
  select(timepoint, de_thresh, sensitivity, specificity, ppv, fnr) %>%
  pivot_longer(c(sensitivity, specificity, ppv, fnr),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric,
                         levels = c("sensitivity", "specificity", "ppv", "fnr"),
                         labels = c("Sensitivity\n(caught true misses)",
                                    "Specificity\n(spared true hits)",
                                    "PPV\n(flagged = truly missed)",
                                    "FNR\n(missed misses)"))) %>%
  ggplot(aes(x = factor(de_thresh), y = timepoint, fill = value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(value, 2)), size = 3.5, color = "white", fontface = "bold") +
  facet_wrap(~metric, ncol = 2) +
  scale_fill_gradient2(low = "#CC4444", mid = "#E09B2B", high = "#2C6E8A",
                       midpoint = 0.5, limits = c(0, 1), name = NULL) +
  labs(
    title    = paste0("Redosing Decision Rule Metrics  |  CME threshold = ", round(cme_cutoff * 100), "%"),
    subtitle = "Redose if drug effect ≤ threshold. Positive = missed CME.",
    x        = "Drug effect threshold (redose if ≤ this)",
    y        = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))

#------------------------------------------------------------------------------------------------------
# DECISION RULE ROC: 90 MIN DRUG EFFECT → MISSED CME
# Direct ROC — drug effect as score, missed CME as outcome
# Higher drug effect = lower risk of missing CME, so direction = ">"
#------------------------------------------------------------------------------------------------------

df_90 <- df %>%
  filter(timepoint == 90) %>%
  group_by(study_id, ID, session, meq.total.pctmax) %>%
  summarise(drug_effect = mean(Avg_OverallDrugEffect, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(drug_effect)) %>%
  mutate(missed_cme = as.integer(meq.total.pctmax < cme_cutoff))

roc_90  <- roc(df_90$missed_cme, df_90$drug_effect, direction = ">", quiet = TRUE)
auc_90  <- round(as.numeric(auc(roc_90)), 3)
ci_90   <- ci.auc(roc_90, conf.level = 0.95)
op      <- coords(roc_90, x = 1, input = "threshold",
                  ret = c("sensitivity", "specificity"), transpose = FALSE)

par(mar = c(4, 4, 4, 1))
plot(roc_90,
     col  = "#2C6E8A", lwd = 2.5,
     xlim = c(1, 0), ylim = c(0, 1),
     main = paste0("ROC: 90 min Drug Effect → Missed CME  |  MEQ < ", round(cme_cutoff * 100), "%\n",
                   "AUC = ", auc_90, "  (95% CI: ", round(ci_90[1], 3), "–", round(ci_90[3], 3), ")"))
abline(a = 1, b = -1, col = "#DDDDDD", lty = 2)
points(op$specificity, op$sensitivity, pch = 19, cex = 2, col = "#E05A2B")
text(op$specificity - 0.05, op$sensitivity - 0.05,
     labels = paste0("Threshold ≤ 1\nSens=", round(op$sensitivity, 2),
                     ", Spec=", round(op$specificity, 2)),
     cex = 0.85, col = "#E05A2B", adj = c(1, 1))

cat("\nDecision rule: redose if 90 min drug effect ≤ 1\n")
cat("AUC:         ", auc_90, "\n")
cat("95% CI:      ", round(ci_90[1], 3), "–", round(ci_90[3], 3), "\n")
cat("Sensitivity: ", round(op$sensitivity, 3), "\n")
cat("Specificity: ", round(op$specificity, 3), "\n")
cat("FPR:         ", round(1 - op$specificity, 3), "\n")
cat("FNR:         ", round(1 - op$sensitivity, 3), "\n")

#------------------------------------------------------------------------------------------------------
# ROBUSTNESS: DECISION RULE ACROSS MEQ THRESHOLDS
# Shows that 90 min ≤ 1 performs consistently regardless of where CME line is drawn
# Loops over MEQ cutoffs (40–80%) x drug effect thresholds (integer only, 0–3) at 90 min only
#------------------------------------------------------------------------------------------------------

meq_cutoffs <- seq(0.40, 0.80, by = 0.10)
de_thresh_int <- c(0, 1, 2, 3)

robust_table <- map_dfr(meq_cutoffs, function(meq_cut) {
  map_dfr(de_thresh_int, function(thresh) {
    
    d <- df %>%
      filter(timepoint == 90) %>%
      group_by(study_id, ID, session, meq.total.pctmax) %>%
      summarise(drug_effect = mean(Avg_OverallDrugEffect, na.rm = TRUE), .groups = "drop") %>%
      filter(!is.na(drug_effect)) %>%
      mutate(
        missed_cme = as.integer(meq.total.pctmax < meq_cut),
        flagged    = as.integer(drug_effect <= thresh)
      )
    
    tp_count <- sum(d$missed_cme == 1 & d$flagged == 1)
    tn_count <- sum(d$missed_cme == 0 & d$flagged == 0)
    fp_count <- sum(d$missed_cme == 0 & d$flagged == 1)
    fn_count <- sum(d$missed_cme == 1 & d$flagged == 0)
    n_missed <- sum(d$missed_cme)
    n_hit    <- nrow(d) - n_missed
    
    tibble(
      meq_cutoff   = meq_cut,
      de_thresh    = thresh,
      n_missed_cme = n_missed,
      pct_missed   = round(n_missed / nrow(d), 2),
      sensitivity  = round(tp_count / n_missed,              3),
      specificity  = round(tn_count / n_hit,                 3),
      ppv          = round(tp_count / (tp_count + fp_count), 3),
      fpr          = round(fp_count / n_hit,                 3),
      fnr          = round(fn_count / n_missed,              3)
    )
  })
}) %>%
  mutate(across(c(sensitivity, specificity, ppv, fpr, fnr),
                ~ ifelse(is.nan(.), NA_real_, .)))

print(robust_table, n = Inf)

# Faceted heatmap: one facet per drug effect threshold, FPR and sensitivity across MEQ cutoffs
robust_table %>%
  select(meq_cutoff, de_thresh, sensitivity, fpr) %>%
  pivot_longer(c(sensitivity, fpr), names_to = "metric", values_to = "value") %>%
  mutate(
    metric    = factor(metric,
                       levels = c("sensitivity", "fpr"),
                       labels = c("Sensitivity\n(caught true misses)",
                                  "FPR\n(unnecessary redoses)")),
    de_thresh = factor(paste0("Drug effect ≤ ", de_thresh))
  ) %>%
  ggplot(aes(x = factor(round(meq_cutoff * 100)), y = metric, fill = value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(value, 2)), size = 3.8, color = "white", fontface = "bold") +
  facet_wrap(~de_thresh, ncol = 2) +
  scale_fill_gradient2(low = "#CC4444", mid = "#E09B2B", high = "#2C6E8A",
                       midpoint = 0.5, limits = c(0, 1), name = NULL) +
  labs(
    title    = "Decision Rule Robustness Across MEQ Thresholds  |  90 min",
    subtitle = "Each panel = one drug effect threshold. Columns = MEQ cutoff definition.",
    x        = "MEQ threshold (% max)",
    y        = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))
  




#------------------------------------------------------------------------------------------------------
# SPEARMAN CORRELATION: DRUG EFFECT VS MEQ BY TIMEPOINT
#------------------------------------------------------------------------------------------------------

cor_table <- map_dfr(c(30, 60, 90, 120, 180, 240, 360), function(tp) {
  d <- df %>%
    filter(timepoint == tp) %>%
    group_by(study_id, ID, session, meq.total.pctmax) %>%
    summarise(drug_effect = mean(Avg_OverallDrugEffect, na.rm = TRUE), .groups = "drop") %>%
    filter(!is.na(drug_effect), !is.na(meq.total.pctmax))
  
  ct <- cor.test(d$drug_effect, d$meq.total.pctmax, method = "spearman", exact = FALSE)
  
  tibble(
    timepoint = paste0(tp, " min"),
    n         = nrow(d),
    rho       = round(ct$estimate, 3),
    p         = round(ct$p.value, 4),
    p_fmt     = ifelse(ct$p.value < 0.001, "<0.001", as.character(round(ct$p.value, 3)))
  )
})

print(cor_table)

# Panel scatterplot
df %>%
  group_by(study_id, ID, session, timepoint, meq.total.pctmax) %>%
  summarise(drug_effect = mean(Avg_OverallDrugEffect, na.rm = TRUE), .groups = "drop") %>%
  left_join(cor_table %>%
              mutate(timepoint = as.numeric(str_extract(timepoint, "\\d+"))),
            by = "timepoint") %>%
  mutate(tp_label  = paste0(timepoint, " min"),
         rho_label = paste0("rho = ", rho, ", p ", p_fmt)) %>%
  ggplot(aes(x = drug_effect, y = meq.total.pctmax)) +
  geom_jitter(alpha = 0.2, width = 0.05, size = 1.5, color = "#2C6E8A") +
  geom_smooth(method = "lm", se = TRUE, color = "#E05A2B",
              fill = "#E05A2B", alpha = 0.15) +
  geom_text(aes(label = rho_label), x = 0.2, y = 0.95,
            hjust = 0, size = 3.8, color = "#333333", fontface = "italic",
            check_overlap = TRUE) +
  facet_wrap(~tp_label, nrow = 1) +
  scale_x_continuous(breaks = 0:4) +
  scale_y_continuous(limits = c(-0.05, 1.05), labels = scales::percent) +
  labs(
    title = "Drug Effect vs MEQ Total % Max by Timepoint",
    x     = "Avg drug effect rating (0-4)",
    y     = "MEQ total % max"
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        strip.text       = element_text(face = "bold"))


redose_table %>%
  filter(timepoint == "90 min") %>%
  select(de_thresh, ppv, fpr) %>%
  pivot_longer(c(ppv, fpr), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric,
                         levels = c("ppv", "fpr"),
                         labels = c("PPV (redose was correct)",
                                    "FPR (unnecessary redose)"))) %>%
  ggplot(aes(x = de_thresh, y = value, color = metric)) +
  geom_line(linewidth = 1.4) +
  geom_point(size = 3.5) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#555555") +
  annotate("text", x = 1.08, y = 0.95, label = "Chosen rule (≤1)",
           hjust = 0, size = 3.8, color = "#555555") +
  scale_color_manual(values = c("PPV (redose was correct)"   = "#2C6E8A",
                                "FPR (unnecessary redose)" = "#E05A2B"),
                     name = NULL) +
  scale_x_continuous(breaks = seq(0, 3, by = 0.5)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(
    title    = "PPV and FPR by Drug Effect Threshold  |  90 min",
    subtitle = paste0("CME threshold = ", round(cme_cutoff * 100), "%  |  Redose if drug effect ≤ threshold"),
    x        = "Drug effect threshold (redose if ≤ this)",
    y        = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "top")


map_dfr(seq(0.40, 1.00, by = 0.10), function(cut) {
  df %>%
    filter(timepoint == 90) %>%
    group_by(study_id, ID, session, meq.total.pctmax) %>%
    summarise(drug_effect = mean(Avg_OverallDrugEffect, na.rm = TRUE),
              .groups = "drop") %>%
    summarise(
      meq_cutoff   = cut,
      n_total      = n(),
      n_cme        = sum(meq.total.pctmax >= cut),
      pct_cme      = round(mean(meq.total.pctmax >= cut) * 100, 1),
      mean_meq     = round(mean(meq.total.pctmax), 3),
      mean_de_hit  = round(mean(drug_effect[meq.total.pctmax >= cut],  na.rm = TRUE), 2),
      mean_de_miss = round(mean(drug_effect[meq.total.pctmax < cut],   na.rm = TRUE), 2)
    )
}) %>%
  print(n = Inf)