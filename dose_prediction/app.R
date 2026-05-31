## PPV / FPR Decision Rule Explorer
## nsepeda1@jh.edu

library(here)
library(tidyverse)
library(shiny)
library(bslib)

#------------------------------------------------------------------------------------------------------
# LOAD DATA
#------------------------------------------------------------------------------------------------------

monrat_avg <- read_csv(here("data", "Combined_MonitorSessionRating_AVG.csv"),
                       show_col_types = FALSE)
guidefx    <- read_csv(here("data", "GuideEffectsData_CLEAN.csv"),
                       show_col_types = FALSE)

guidefx_clean <- guidefx %>%
  select(study_id, study_name, ID, session, meq.total.pctmax) %>%
  mutate(across(c(study_id, ID, session), as.numeric))

df <- monrat_avg %>%
  filter(timepoint %in% c(30, 60, 90, 120)) %>%
  select(study_id, study_name, ID, session, timepoint, Avg_OverallDrugEffect) %>%
  inner_join(guidefx_clean, by = c("study_id", "study_name", "ID", "session")) %>%
  filter(!is.na(meq.total.pctmax), !is.na(Avg_OverallDrugEffect)) %>%
  mutate(across(c(study_id, ID, session, timepoint), as.numeric))

#------------------------------------------------------------------------------------------------------
# UI
#------------------------------------------------------------------------------------------------------

ui <- page_sidebar(
  title = "Redosing Decision Rule Explorer",
  theme = bs_theme(
    bootswatch = "flatly",
    base_font  = font_google("IBM Plex Sans"),
    code_font  = font_google("IBM Plex Mono"),
    primary    = "#2C6E8A",
    bg         = "#F7F9FB",
    fg         = "#1A1A2E"
  ),
  
  sidebar = sidebar(
    width = 260,
    
    tags$h6("TIMEPOINT", style = "letter-spacing:0.1em; color:#888; margin-bottom:4px;"),
    radioButtons("timepoint", label = NULL,
                 choices  = c("30 min" = 30, "60 min" = 60,
                              "90 min" = 90, "120 min" = 120),
                 selected = 90),
    
    hr(),
    
    tags$h6("MEQ THRESHOLD", style = "letter-spacing:0.1em; color:#888; margin-bottom:4px;"),
    sliderInput("meq_cutoff", label = NULL,
                min = 0.30, max = 0.85, value = 0.60, step = 0.05,
                post = " (% max)"),
    
    hr(),
    
    tags$h6("CHOSEN DRUG EFFECT CUTOFF",
            style = "letter-spacing:0.1em; color:#888; margin-bottom:4px;"),
    sliderInput("chosen_thresh", label = NULL,
                min = 0, max = 3, value = 1, step = 0.5),
    tags$small(style = "color:#888;",
               "Vertical line — redose if drug effect ≤ this value")
  ),
  
  layout_columns(
    col_widths = c(8, 4),
    
    card(
      card_header("PPV & FPR by Drug Effect Threshold"),
      plotOutput("ppv_fpr_plot", height = "420px")
    ),
    
    card(
      card_header("At Chosen Cutoff"),
      tableOutput("metrics_table"),
      hr(),
      uiOutput("summary_badges"),
      hr(),
      tags$div(
        style = "padding: 8px 4px; font-size: 0.85em; color: #555;",
        tags$p(style = "font-weight: 600; margin-bottom: 6px;",
               "Decision Rule Goals"),
        tags$p(style = "margin-bottom: 6px;",
               tags$span(style = "color: #E05A2B; font-weight:600;", "FPR ≤ 10% "),
               "— fewer than 1 in 10 redoses should go to someone who would have",
               " reached MEQ goal without intervention."),
        tags$p(style = "margin-bottom: 0;",
               tags$span(style = "color: #2C6E8A; font-weight:600;", "PPV ≥ 80% "),
               "— at least 8 in 10 redose decisions should be correct.",
               " When we intervene, we should be confident the person needed it.")
      )
    )
  )
)

#------------------------------------------------------------------------------------------------------
# SERVER
#------------------------------------------------------------------------------------------------------

server <- function(input, output, session) {
  
  # Compute redose table for selected timepoint + MEQ cutoff
  redose_data <- reactive({
    tp         <- as.numeric(input$timepoint)
    meq_cutoff <- input$meq_cutoff
    
    map_dfr(seq(0, 3, by = 0.5), function(thresh) {
      d <- df %>%
        filter(timepoint == tp) %>%
        group_by(study_id, ID, session, meq.total.pctmax) %>%
        summarise(drug_effect = mean(Avg_OverallDrugEffect, na.rm = TRUE),
                  .groups = "drop") %>%
        filter(!is.na(drug_effect)) %>%
        mutate(
          missed_cme = as.integer(meq.total.pctmax < meq_cutoff),
          flagged    = as.integer(drug_effect <= thresh)
        )
      
      tp_count <- sum(d$missed_cme == 1 & d$flagged == 1)
      tn_count <- sum(d$missed_cme == 0 & d$flagged == 0)
      fp_count <- sum(d$missed_cme == 0 & d$flagged == 1)
      fn_count <- sum(d$missed_cme == 1 & d$flagged == 0)
      n_missed <- sum(d$missed_cme)
      n_hit    <- nrow(d) - n_missed
      
      tibble(
        de_thresh   = thresh,
        n_sessions  = nrow(d),
        n_missed    = n_missed,
        n_flagged   = tp_count + fp_count,
        sensitivity = tp_count / n_missed,
        specificity = tn_count / n_hit,
        ppv         = tp_count / (tp_count + fp_count),
        npv         = tn_count / (tn_count + fn_count),
        fpr         = fp_count / n_hit,
        fnr         = fn_count / n_missed
      )
    }) %>%
      mutate(across(c(sensitivity, specificity, ppv, npv, fpr, fnr),
                    ~ ifelse(is.nan(.), NA_real_, .)))
  })
  
  # PPV / FPR plot
  output$ppv_fpr_plot <- renderPlot({
    d      <- redose_data()
    thresh <- input$chosen_thresh
    meq_pct <- round(input$meq_cutoff * 100)
    tp      <- input$timepoint
    
    d_long <- d %>%
      select(de_thresh, ppv, fpr) %>%
      pivot_longer(c(ppv, fpr), names_to = "metric", values_to = "value") %>%
      mutate(metric = factor(metric,
                             levels = c("ppv", "fpr"),
                             labels = c("PPV (redose was correct)",
                                        "FPR (unnecessary redose)")))
    
    # Values at chosen threshold
    at_thresh <- d %>% filter(de_thresh == thresh)
    
    ggplot(d_long, aes(x = de_thresh, y = value, color = metric)) +
      geom_line(linewidth = 1.6) +
      geom_point(size = 4) +
      geom_vline(xintercept = thresh, linetype = "dashed",
                 color = "#444444", linewidth = 0.8) +
      annotate("text", x = thresh + 0.07, y = 0.97,
               label = paste0("Cutoff ≤ ", thresh),
               hjust = 0, size = 4, color = "#444444") +
      # Annotate PPV and FPR values at chosen threshold
      {if (nrow(at_thresh) > 0) list(
        annotate("point", x = thresh, y = at_thresh$ppv,
                 size = 6, shape = 21, fill = "#2C6E8A", color = "white"),
        annotate("point", x = thresh, y = at_thresh$fpr,
                 size = 6, shape = 21, fill = "#E05A2B", color = "white")
      )} +
      scale_color_manual(
        values = c("PPV (redose was correct)" = "#2C6E8A",
                   "FPR (unnecessary redose)" = "#E05A2B"),
        name = NULL
      ) +
      scale_x_continuous(breaks = seq(0, 3, by = 0.5)) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
      labs(
        x        = "Drug effect threshold (redose if ≤ this)",
        y        = NULL,
        subtitle = paste0(tp, " min  |  Threshold = MEQ ≥ ", meq_pct, "%  |  N = ",
                          d$n_sessions[1], " sessions  (",
                          d$n_missed[1], " missed threshold)")
      ) +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position  = "top",
        plot.background  = element_rect(fill = "#F7F9FB", color = NA),
        plot.subtitle    = element_text(color = "#888888", size = 11)
      )
  })
  
  # Metrics table at chosen threshold
  output$metrics_table <- renderTable({
    d      <- redose_data()
    thresh <- input$chosen_thresh
    row    <- d %>% filter(de_thresh == thresh)
    req(nrow(row) > 0)
    
    tibble(
      Metric = c("Drug effect threshold", "N sessions", "N missed threshold",
                 "N flagged for redose",
                 "PPV", "FPR", "Sensitivity", "Specificity", "NPV", "FNR"),
      Value  = c(
        paste0("≤ ", thresh),
        row$n_sessions,
        row$n_missed,
        row$n_flagged,
        paste0(round(row$ppv * 100, 1), "%"),
        paste0(round(row$fpr * 100, 1), "%"),
        paste0(round(row$sensitivity * 100, 1), "%"),
        paste0(round(row$specificity * 100, 1), "%"),
        paste0(round(row$npv * 100, 1), "%"),
        paste0(round(row$fnr * 100, 1), "%")
      )
    )
  }, striped = TRUE, hover = TRUE, bordered = FALSE, spacing = "s")
  
  # Summary badges
  output$summary_badges <- renderUI({
    d      <- redose_data()
    thresh <- input$chosen_thresh
    row    <- d %>% filter(de_thresh == thresh)
    req(nrow(row) > 0)
    
    ppv_col <- ifelse(row$ppv >= 0.85, "#2C6E8A",
                      ifelse(row$ppv >= 0.70, "#E09B2B", "#CC4444"))
    fpr_col <- ifelse(row$fpr <= 0.10, "#2C6E8A",
                      ifelse(row$fpr <= 0.20, "#E09B2B", "#CC4444"))
    
    tags$div(
      style = "display:flex; gap:12px; justify-content:center; padding:8px 0;",
      tags$div(
        style = "text-align:center;",
        tags$div(style = paste0("font-size:2.2em; font-weight:700; color:", ppv_col),
                 paste0(round(row$ppv * 100), "%")),
        tags$div(style = "color:#888; font-size:0.8em;", "PPV")
      ),
      tags$div(
        style = "text-align:center;",
        tags$div(style = paste0("font-size:2.2em; font-weight:700; color:", fpr_col),
                 paste0(round(row$fpr * 100), "%")),
        tags$div(style = "color:#888; font-size:0.8em;", "FPR")
      )
    )
  })
}

#------------------------------------------------------------------------------------------------------
# RUN
#------------------------------------------------------------------------------------------------------

shinyApp(ui, server)