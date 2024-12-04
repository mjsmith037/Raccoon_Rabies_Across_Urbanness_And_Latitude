library(GLMMadaptive) # for statistical models
library(splines)      # for statistical models
library(pROC)         # for model fit
library(kableExtra)   # for tables
library(ggbeeswarm)   # for figures
library(ggpubr)       # for figures
library(patchwork)    # for figures
library(magrittr)     # for data wrangling
library(tidyverse)    # for data wrangling

#### settings and ancillary functions ####
theme_set(theme_bw())
base_tests <- 2^0
num_samples <- 100
resolution <- 100
monthly_data <- read_csv("../data/hypothetical_data.csv")
persistence_data <- read_csv("../data/persistence_5.11_data.csv")
my_colours <- c("#00A676", "#934D77", "#485696", "#FC7A1E")

get_sig <- function(p) {
  case_when(p < 0.001 ~ "***",
            p < 0.01 ~ "**",
            p < 0.05 ~ "*",
            TRUE ~ "")
}

#### best fitting models ####
positivity_model_name <- "density_1_1_month_6_1_latitude_5_1"
load(str_glue("../results/positivity_{positivity_model_name}.RData"))
positivity_model_name <- str_c("positivity_", positivity_model_name)
positivity_model <- get(positivity_model_name)

presence_model_name <- "density_1_3_month_7_1_latitude_7_1"
load(str_glue("../results/presence_{presence_model_name}.RData"))
presence_model_name <- str_c("presence_", presence_model_name)
presence_model <- get(presence_model_name)

n <- 5; m <- 11
persistence_model_name <- "density_1_1_month_1_3_latitude_1_2_5_11"
load(str_glue("../results/persistence_{persistence_model_name}.RData"))
persistence_model_name <- str_glue("persistence_{n}.{m}_{str_remove(persistence_model_name, str_glue('_{n}_{m}$'))}")
persistence_model <- get(persistence_model_name)

#### predictions ####
predictions_month <- tibble(
  density      = sample(monthly_data$density, num_samples, replace=TRUE),
  spatial_lag  = sample(monthly_data$spatial_lag, num_samples, replace=TRUE),
  temporal_lag = sample(monthly_data$temporal_lag, num_samples, replace=TRUE),
  habitat      = sample(monthly_data$habitat, num_samples, replace=TRUE),
  total_tested = sample(monthly_data$total_tested[monthly_data$total_tested != 0], num_samples, replace=TRUE),
  total_tested_over_lag = sample(persistence_data$total_tested_over_lag[persistence_data$total_tested_over_lag != 0],
                                 num_samples, replace=TRUE)) %>%
  crossing(month  = seq(min(monthly_data$month), max(monthly_data$month), length.out=resolution),
           latitude   = quantile(monthly_data$latitude, probs=c(0.11, 0.55, 0.83))) %>%
  group_split(month, .keep=TRUE) %>%
  lapply(function(dat) {
    full_join(effectPlotData(positivity_model, newdata=dat),
              effectPlotData(presence_model, newdata=dat),
              by=c("month", "density", "latitude",
                   "spatial_lag", "temporal_lag",
                   "habitat", "total_tested", "total_tested_over_lag"),
              suffix=c("", "_zi")) %>%
      full_join(effectPlotData(persistence_model, newdata=dat),
                by=c("month", "density", "latitude",
                     "spatial_lag", "temporal_lag",
                     "habitat", "total_tested", "total_tested_over_lag"),
                suffix=c("", "_end"))}) %>%
  bind_rows() %>%
  # un-transform variables after getting predictions
  mutate(across(c(pred, low, upp, pred_zi, low_zi, upp_zi, pred_end, low_end, upp_end), plogis))

predictions_density <- tibble(
  month        = sample(monthly_data$month, num_samples, replace=TRUE),
  spatial_lag  = sample(monthly_data$spatial_lag, num_samples, replace=TRUE),
  temporal_lag = sample(monthly_data$temporal_lag, num_samples, replace=TRUE),
  habitat      = sample(monthly_data$habitat, num_samples, replace=TRUE),
  total_tested = sample(monthly_data$total_tested[monthly_data$total_tested != 0], num_samples, replace=TRUE),
  total_tested_over_lag = sample(persistence_data$total_tested_over_lag[persistence_data$total_tested_over_lag != 0],
                                 num_samples, replace=TRUE)) %>%
  crossing(density    = seq(min(monthly_data$density), max(monthly_data$density), length.out=resolution),
           latitude   = quantile(monthly_data$latitude, probs=c(0.11, 0.55, 0.83))) %>%
  group_split(density, .keep=TRUE) %>%
  lapply(function(dat) {
    full_join(effectPlotData(positivity_model, newdata=dat),
              effectPlotData(presence_model, newdata=dat),
              by=c("month", "density", "latitude",
                   "spatial_lag", "temporal_lag",
                   "habitat", "total_tested", "total_tested_over_lag"),
              suffix=c("", "_zi")) %>%
      full_join(effectPlotData(persistence_model, newdata=dat),
                by=c("month", "density", "latitude",
                     "spatial_lag", "temporal_lag",
                     "habitat", "total_tested", "total_tested_over_lag"),
                suffix=c("", "_end"))}) %>%
  bind_rows() %>%
  # un-transform variables after getting predictions
  mutate(across(c(pred, low, upp, pred_zi, low_zi, upp_zi, pred_end, low_end, upp_end), plogis))

#### SUMMARY TABLE FOR BEST PREVALENCE MODEL  # # # # # # # # # # # # # # # ####
summary(positivity_model)$coef_table %>%
  as_tibble(rownames="Term") %>%
  select(-`z-value`) %>%
  mutate(sig = get_sig(`p-value`),
         Term = str_remove_all(Term, "\\d$|\\d(?=:)"),
         Term = str_remove_all(Term, "\\(month\\)"),
         Term = str_remove_all(Term, ", df = 1, degree = 1\\)"),
         Term = str_remove_all(Term, "bs\\((?=\\w+($|:))"),
         Term = str_replace_all(Term, c("total_tested"="Total Number of Samples Submitted",
                                        "temporal_lag"="Temporal Lag",
                                        "spatial_lag"="Spatial Effect",
                                        "latitude"="Latitude",
                                        "density"="Population Density",
                                        "habitat"="Raccoon Favorable Habitat",
                                        "month"="Month"))) %>%
  rename(`\\emph{p} value` = `p-value`,
         `Standard Error` = Std.Err,
         ` ` = sig) %>%
  kable(booktabs=TRUE, digits=4, format="latex", escape=FALSE, align="lrrrc", linesep="") %>%
  add_header_above(c(" "=1, "Prevalence"=4)) %>%
  collapse_rows(1, latex_hline="major") %T>%
  save_kable(file="../figures/best_prevalence_model.png")

#### SUMMARY TABLE FOR BEST PRESENCE MODEL # # # # # # # # # # # # # # # # ####
summary(presence_model)$coef_table %>%
  as_tibble(rownames="Term") %>%
  select(-`z-value`) %>%
  mutate(sig = get_sig(`p-value`),
         Term = str_remove_all(Term, "\\d$|\\d(?=:)"),
         Term = str_remove_all(Term, "\\(month\\)"),
         Term = str_remove_all(Term, ", df = 1, degree = 1\\)"),
         Term = str_remove_all(Term, "bs\\((?=\\w+($|:))"),
         Term = str_replace_all(Term, c("total_tested"="Total Number of Samples Submitted",
                                        "temporal_lag"="Temporal Lag",
                                        "spatial_lag"="Spatial Effect",
                                        "latitude"="Latitude",
                                        "density"="Population Density",
                                        "habitat"="Raccoon Favorable Habitat",
                                        "month"="Month"))) %>%
  rename(`\\emph{p} value` = `p-value`,
         `Standard Error` = Std.Err,
         ` ` = sig) %>%
  kable(booktabs=TRUE, digits=4, format="latex", escape=FALSE, align="lrrrc", linesep="") %>%
  add_header_above(c(" "=1, "Presence"=4)) %>%
  collapse_rows(1, latex_hline="major") %T>%
  save_kable(file="../figures/best_presence_model.png")

#### SUMMARY TABLE FOR BEST PERSISTENCE MODEL # # # # # # # # # # # # # # # ####
summary(persistence_model)$coef_table %>%
  as_tibble(rownames="Term") %>%
  select(-`z-value`) %>%
  mutate(sig = get_sig(`p-value`),
         Term = str_remove_all(Term, "\\d$|\\d(?=:)"),
         Term = str_remove_all(Term, "\\(month\\)"),
         Term = str_remove_all(Term, ", df = 1, degree = 1\\)"),
         Term = str_remove_all(Term, "bs\\((?=\\w+($|:))"),
         Term = str_replace_all(Term, c("total_tested_over_lag"="Total Number of Samples Submitted Over Previous 11 Months",
                                        "temporal_lag"="Temporal Lag",
                                        "spatial_lag"="Spatial Effect",
                                        "latitude"="Latitude",
                                        "density"="Population Density",
                                        "habitat"="Raccoon Favorable Habitat",
                                        "month"="Month"))) %>%
  rename(`\\emph{p} value` = `p-value`,
         `Standard Error` = Std.Err,
         ` ` = sig) %>%
  kable(booktabs=TRUE, digits=4, format="latex", escape=FALSE, align="lrrrc", linesep="") %>%
  add_header_above(c(" "=1, "Persistence"=4)) %>%
  collapse_rows(1, latex_hline="major") %T>%
  save_kable(file="../figures/best_persistence_model.png")

#### TESTS AND POSITIVITY BY MONTH  # # # # # # # # # # # # # # # # # # # # ####
monthly_data %>%
  mutate(region = case_when(latitude <= quantile(unique(latitude), 0.25) ~
                              "South (lower quartile by latitude)",
                            latitude >= quantile(unique(latitude), 0.75) ~
                              "North (upper quartile by latitude)",
                            TRUE ~ "Middle") %>%
           factor(levels=c("South (lower quartile by latitude)",
                           # "Middle",
                           "North (upper quartile by latitude)"))) %>%
  na.omit() %>%
  mutate(positivity = positive / total_tested) %>%
  {ggplot(.) +
      aes(x=month, y=total_tested) +
      geom_quasirandom(aes(group=month), varwidth=TRUE, size=0.25) +
      geom_smooth(colour=my_colours[4], method="gam", formula=y~s(x, bs="cs")) +
      scale_y_log10(breaks=c(0:5, 10*1:5, 100)) +
      ylab("Samples submitted") +
      facet_wrap(~region, nrow=1) +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      ggplot(.) +
      aes(x=month, y=positivity) +
      geom_quasirandom(aes(group=month), varwidth=TRUE, size=0.25) +
      geom_smooth(colour=my_colours[4], method="gam", formula=y~s(x, bs="cs")) +
      scale_x_continuous(name=NULL, breaks=1:12, labels=month.name) +
      scale_y_continuous(name="Percent Positivity", limits=0:1, breaks=0:4/4,
                         labels=c("0%", "25%", "50%", "75%", "100%")) +
      facet_wrap(~region, nrow=1) +
      theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1),
            strip.text=element_blank(),
            strip.background=element_blank()) +
      plot_layout(ncol=1)}
ggsave(filename="../figures/tests_and_positivity.png", width=8, height=6)

#### MONTH LATITUDE INTERACTION PLOT  # # # # # # # # # # # # # # # # # # # ####
{predictions_month %>%
    group_by(month, latitude) %>%
    summarise(low_zi = quantile(pred_zi, 0.1), high_zi=quantile(pred_zi, 0.9), pred_zi = mean(pred_zi), .groups="drop") %>%
    mutate(latitude = latitude %>% round() %>% str_c("°") %>%
             factor(ordered=TRUE, labels=c("Southern", "Mid-latitude", "Northern"))) %>%
    {ggplot(.) +
        aes(x=month, group=str_c(latitude)) +
        geom_ribbon(aes(ymin=low_zi, ymax=high_zi, fill=latitude), alpha=0.2) +
        geom_smooth(aes(y=pred_zi, colour=latitude), linewidth=1, se=FALSE, span=0.2) +
        coord_cartesian(ylim=c(0,1)) +
        scale_y_continuous(name="Presence", breaks=0:4/4) +
        scale_x_continuous(name=NULL, breaks=1:12, labels=month.name) +
        scale_colour_manual(name="Latitude:", values=my_colours,
                            aesthetics=c("colour", "fill")) +
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank())}} +
  
  {predictions_month %>%
      group_by(month, latitude) %>%
      summarise(low = quantile(pred, 0.1), high=quantile(pred, 0.9), pred = mean(pred), .groups="drop") %>%
      mutate(latitude = latitude %>% round() %>% str_c("°") %>%
               factor(ordered=TRUE, labels=c("Southern", "Mid-latitude", "Northern"))) %>%
      {ggplot(.) +
          aes(x=month) +
          geom_ribbon(aes(ymin=low, ymax=high, fill=latitude), alpha=0.2) +
          geom_smooth(aes(y=pred, colour=latitude), linewidth=1, se=FALSE, span=0.2) +
          coord_cartesian(ylim=c(0,1)) +
          scale_y_continuous(name="Prevalence", breaks=0:4/4) +#,
                             # labels=c("0%", "25%", "50%", "75%", "100%")) +
          scale_x_continuous(name=NULL, breaks=1:12, labels=month.name) +
          scale_colour_manual(name="Latitude:", values=my_colours,
                              aesthetics=c("colour", "fill")) +
          theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1))}} +
  
  {predictions_month %>%
      group_by(month, latitude) %>%
      summarise(low_end = quantile(pred_end, 0.1), high_end = quantile(pred_end, 0.9),
                pred_end = mean(pred_end), .groups="drop") %>%
      mutate(latitude = latitude %>% round() %>% str_c("°") %>%
               factor(ordered=TRUE, labels=c("Southern", "Mid-latitude", "Northern"))) %>%
      {ggplot(.) +
          aes(x=month) +
          geom_ribbon(aes(ymin=low_end, ymax=high_end, fill=latitude), alpha=0.2) +
          geom_smooth(aes(y=pred_end, colour=latitude), linewidth=1, se=FALSE, span=0.2) +
          coord_cartesian(ylim=c(0,1)) +
          scale_y_continuous(name="Persistence", breaks=0:4/4) +
          scale_x_continuous(name=NULL, breaks=1:12, labels=month.name) +
          scale_colour_manual(name="Latitude:", values=my_colours,
                              aesthetics=c("colour", "fill")) +
          theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1))}} +
  patchwork::plot_layout(ncol=1, guides="collect") &
  theme(legend.position="bottom")

ggsave(filename="../figures/month_latitude_interaction.png", width=5, height=6)

#### DENSITY LATITUDE INTERACTION PLOT  # # # # # # # # # # # # # # # # # # ####
{predictions_density %>%
    group_by(density, latitude) %>%
    summarise(low_zi = quantile(pred_zi, 0.1), high_zi=quantile(pred_zi, 0.9), pred_zi = mean(pred_zi), .groups="drop") %>%
    mutate(latitude = latitude %>% round() %>% str_c("°") %>%
             factor(ordered=TRUE, labels=c("Southern", "Mid-latitude", "Northern"))) %>%
    ggplot() +
    aes(x=density) +
    geom_ribbon(aes(ymin=low_zi, ymax=high_zi, fill=latitude), alpha=0.2) +
    geom_smooth(aes(y=pred_zi, colour=latitude), linewidth=1, se=FALSE) +
    coord_cartesian(ylim=c(0,1)) +
    scale_y_continuous(name="Presence", breaks=0:4/4) +
    scale_colour_manual(name="Latitude:", values=my_colours,
                        aesthetics=c("colour", "fill")) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank())} +
  
  {predictions_density %>%
      group_by(density, latitude) %>%
      summarise(low = quantile(pred, 0.1), high=quantile(pred, 0.9), pred = mean(pred), .groups="drop") %>%
      mutate(latitude = latitude %>% round() %>% str_c("°") %>%
               factor(ordered=TRUE, labels=c("Southern", "Mid-latitude", "Northern"))) %>%
      ggplot() +
      aes(x=density) +
      geom_ribbon(aes(ymin=low, ymax=high, fill=latitude), alpha=0.2) +
      geom_smooth(aes(y=pred, colour=latitude), linewidth=1, se=FALSE) +
      coord_cartesian(ylim=c(0,1)) +
      scale_y_continuous(name="Prevalence", breaks=0:4/4) +#,
                         # labels=c("0%", "25%", "50%", "75%", "100%")) +
      scale_colour_manual(name="Latitude:", values=my_colours,
                          aesthetics=c("colour", "fill")) +
      theme(axis.title.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.x=element_blank())} +
  
  {predictions_density %>%
      group_by(density, latitude) %>%
      summarise(low_end = quantile(pred_end, 0.1), high_end=quantile(pred_end, 0.9), pred_end = mean(pred_end), .groups="drop") %>%
      mutate(latitude = latitude %>% round() %>% str_c("°") %>%
               factor(ordered=TRUE, labels=c("Southern", "Mid-latitude", "Northern"))) %>%
      ggplot() +
      aes(x=density) +
      geom_ribbon(aes(ymin=low_end, ymax=high_end, fill=latitude), alpha=0.2) +
      geom_smooth(aes(y=pred_end, colour=latitude), linewidth=1, se=FALSE) +
      geom_rug(linewidth=0.67, alpha=0.2, data=monthly_data %>% distinct(fips, density)) +
      coord_cartesian(ylim=c(0,1)) +
      scale_y_continuous(name="Persistence", breaks=0:4/4) +
      scale_colour_manual(name="Latitude:", values=my_colours,
                          aesthetics=c("colour", "fill"))} +
  
  plot_layout(ncol=1, guides="collect") &
  scale_x_log10(expand=expansion(), name=expression(paste("Human Population Density (persons/", km^2, ")"))) &
  theme(legend.position="bottom")

ggsave(filename="../figures/latitude_and_density.png", width=5, height=6)

#### HABITAT (NOT SIGNIFICANT)  # # # # # # # # # # # # # # # # # # # # # # ####
cutoff_presence <- 0.16; cutoff_positivity <- 0.17; cutoff_persistence <- 0.16;
(monthly_data %>% mutate(presence=as.integer((positive / total_tested) > 0)) %>%
   {ggplot(.) +
       aes(x=habitat, y=presence) +
       geom_quasirandom(orientation="y", size=0.25, varwidth=TRUE, width=0.1, alpha=0.05) +
       geom_smooth(linewidth=1, colour=my_colours[4], formula=y~x,
                   method = "glm", method.args = list(family = "binomial"),
                   data=filter(., habitat <= cutoff_presence)) +
       stat_regline_equation(data=filter(., habitat <= cutoff_presence),
                             aes(label=paste(..eq.label.., ..adj.rr.label.., sep="~~~~")),
                             colour=my_colours[4], label.x=0.325, label.y=0.85) +
       geom_smooth(linewidth=1, colour=my_colours[3], formula=y~x,
                   method = "glm", method.args = list(family = "binomial"),
                   data=filter(., habitat > cutoff_presence)) +
       stat_regline_equation(data=filter(., habitat > cutoff_presence),
                             aes(label=paste(..eq.label.., ..adj.rr.label.., sep="~~~~")),
                             colour=my_colours[3], label.x=0.325, label.y=0.675) +
       scale_y_continuous(name="Presence", breaks=0:1) +
       theme(axis.text.x=element_blank(),
             axis.title.x=element_blank(),
             axis.ticks.x=element_blank())}) +

  (monthly_data %>%
     mutate(positivity = positive / total_tested) %>%
     {ggplot(.) +
         aes(x=habitat, y=positivity) +
         geom_point(size=0.25, alpha=0.05) +
         geom_smooth(linewidth=1, colour=my_colours[4], formula=y~x,
                     method="lm", data=filter(., habitat <= cutoff_positivity)) +
         stat_regline_equation(data=filter(., habitat <= cutoff_positivity),
                               aes(label=paste(..eq.label.., ..adj.rr.label.., sep="~~~~")),
                               colour=my_colours[4], label.x=0.325, label.y=0.85) +
         geom_smooth(linewidth=1, colour=my_colours[3], formula=y~x,
                     method="lm", data=filter(., habitat > cutoff_positivity)) +
         stat_regline_equation(data=filter(., habitat > cutoff_positivity),
                               aes(label=paste(..eq.label.., ..adj.rr.label.., sep="~~~~")),
                               colour=my_colours[3], label.x=0.325, label.y=0.675) +
         scale_y_continuous(name="Test Positivity", breaks=0:4/4) +#,
                            # labels=c("0%", "25%", "50%", "75%", "100%")) +
         theme(axis.text.x=element_blank(),
               axis.title.x=element_blank(),
               axis.ticks.x=element_blank())}) +

  (persistence_data %>%
     mutate(persistence = as.integer(persistence)) %>%
     {ggplot(.) +
         aes(x=habitat, y=persistence) +
         geom_quasirandom(orientation="y", size=0.25, varwidth=TRUE, width=0.1, alpha=0.05) +
         geom_smooth(linewidth=1, colour=my_colours[4], formula=y~x,
                     method = "glm", method.args = list(family = "binomial"),
                     data=filter(., habitat <= cutoff_persistence)) +
         stat_regline_equation(data=filter(., habitat <= cutoff_persistence),
                               aes(label=paste(after_stat(eq.label), after_stat(adj.rr.label), sep="~~~~")),
                               colour=my_colours[4], label.x=0.325, label.y=0.85) +
         geom_smooth(linewidth=1, colour=my_colours[3], formula=y~x,
                     method = "glm", method.args = list(family = "binomial"),
                     data=filter(., habitat > cutoff_persistence)) +
         stat_regline_equation(data=filter(., habitat > cutoff_persistence),
                               aes(label=paste(after_stat(eq.label), after_stat(adj.rr.label), sep="~~~~")),
                               colour=my_colours[3], label.x=0.325, label.y=0.675) +
         scale_y_continuous(name="Persistence", breaks=0:1)}) +

  plot_layout(ncol=1) &
  scale_x_continuous(name="Proportion of Land-cover Favorable to Raccoons", limits=0:1)
ggsave(filename="../figures/habitat_piecewise_regression.png", width=5, height=7)
