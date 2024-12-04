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
base_tests <- 2^1
monthly_data <- read_csv("../data/hypothetical_data.csv")

get_sig <- function(p) {
  case_when(p < 0.001 ~ "***",
            p < 0.01 ~ "**",
            p < 0.05 ~ "*",
            TRUE ~ "")
}

rug_data <- monthly_data %>%
  select(!c(positive, total_tested)) %>%
  distinct()

#### best fitting models ####
positivity_model_name <- "density_1_1_month_6_1_latitude_5_1"
load(str_glue("../results/{positivity_model_name}.RData"))
positivity_model_name <- str_c("positivity_", positivity_model_name)
positivity_model <- get(positivity_model_name)

detection_model_name <- "density_1_3_month_7_1_latitude_7_1"
load(str_glue("../results/{detection_model_name}.RData"))
detection_model_name <- str_c("detection_", detection_model_name)
detection_model <- get(detection_model_name)

n <- 5; m <- 11
persistence_model_name <- "density_1_1_month_1_3_latitude_1_2_5_11"
load(str_glue("../results/{persistence_model_name}.RData"))
persistence_model_name <- str_glue("endemic_{n}.{m}_{str_remove(persistence_model_name, str_glue('_{n}_{m}$'))}")
persistence_model <- get(persistence_model_name)

#### predictions ####
predictions <- monthly_data %>%
  mutate(latitude = scale(latitude), density = scale(density)) %>%
  {crossing(month        = month.name %>% factor(levels=month.name, ordered=TRUE) %>% as.integer(),
            density      = seq(min(.$density), max(.$density), length.out=10),
            latitude     = seq(min(.$latitude), max(.$latitude), length.out=10),
            spatial_lag  = seq(min(.$spatial_lag), max(.$spatial_lag), length.out=10),
            temporal_lag = seq(min(.$temporal_lag), max(.$temporal_lag), length.out=10),
            habitat      = seq(min(.$habitat), max(.$habitat), length.out=10),
            total_tested = base_tests)} %>%
  group_split(month, .keep=TRUE) %>%
  lapply(. %>% mutate(pred = predict(positivity_model, newdata=.),
                      pred_zi = predict(detection_model, newdata=.))) %>%
  bind_rows() %>%
  # un-transform variables after getting predictions
  mutate(density = density %>% multiply_by(sd(monthly_data$density)) %>% add(mean(monthly_data$density)),
         latitude = latitude %>% multiply_by(sd(monthly_data$latitude)) %>% add(mean(monthly_data$latitude)))

#### SUMMARY TABLE FOR BEST PREVALENCE MODEL  # # # # # # # # # # # # # # # ####
summary(positivity_model)$coef_table %>%
  as_tibble(rownames="Term") %>%
  select(-`z-value`) %>%
  mutate(sig = get_sig(`p-value`),
         Term = str_remove_all(Term, "\\d$|\\d(?=:)"),
         Term = str_remove_all(Term, "\\(month\\)"),
         Term = str_remove_all(Term, ", df = 1, degree = 1\\)"),
         Term = str_remove_all(Term, "bs\\((?=\\w+($|:))"),
         Term = str_replace_all(Term, c("total_tested"="Total Number of Tests",
                                        "temporal_lag"="Temporal Lag",
                                        "spatial_lag"="Spatial Effect",
                                        "latitude"="Latitude",
                                        "density"="Population Density",
                                        "habitat"="Racoon Favorable Habitat",
                                        "month"="Month"))) %>%
  rename(`\\emph{p} value` = `p-value`,
         `Standard Error` = Std.Err,
         ` ` = sig) %>%
  kable(booktabs=TRUE, digits=4, format="latex", escape=FALSE, align="lrrrc") %>%
  add_header_above(c(" "=1, "Percent Positivity"=4)) %>%
  collapse_rows(1, latex_hline="major") %T>%
  save_kable(file="../figures/best_prevalence_model.png")

#### SUMMARY TABLE FOR BEST DETECTION MODEL # # # # # # # # # # # # # # # # ####
summary(detection_model)$coef_table %>%
  as_tibble(rownames="Term") %>%
  select(-`z-value`) %>%
  mutate(sig = get_sig(`p-value`),
         Term = str_remove_all(Term, "\\d$|\\d(?=:)"),
         Term = str_remove_all(Term, "\\(month\\)"),
         Term = str_remove_all(Term, ", df = 1, degree = 1\\)"),
         Term = str_remove_all(Term, "bs\\((?=\\w+($|:))"),
         Term = str_replace_all(Term, c("total_tested"="Total Number of Tests",
                                        "temporal_lag"="Temporal Lag",
                                        "spatial_lag"="Spatial Effect",
                                        "latitude"="Latitude",
                                        "density"="Population Density",
                                        "habitat"="Racoon Favorable Habitat",
                                        "month"="Month"))) %>%
  rename(`\\emph{p} value` = `p-value`,
         `Standard Error` = Std.Err,
         ` ` = sig) %>%
  kable(booktabs=TRUE, digits=4, format="latex", escape=FALSE, align="lrrrc") %>%
  add_header_above(c(" "=1, "Detection"=4)) %>%
  collapse_rows(1, latex_hline="major") %T>%
  save_kable(file="../figures/best_detection_model.png")

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
                                        "habitat"="Racoon Favorable Habitat",
                                        "month"="Month"))) %>%
  rename(`\\emph{p} value` = `p-value`,
         `Standard Error` = `Std.Err`,
         ` ` = sig) %>%
  kable(booktabs=TRUE, digits=4, format="latex", escape=FALSE, align="lrrrc") %>%
  add_header_above(c(" "=1, "Persistence"=4)) %>%
  collapse_rows(1, latex_hline="major") %T>%
  save_kable(file="../figures/best_persistence_model.png")

#### TESTS AND POSITIVITY BY MONTH  # # # # # # # # # # # # # # # # # # # # ####
monthly_data %>%
  mutate(region = case_when(latitude <= quantile(latitude, 0.25) ~
                              "South (lower quartile by latitude)",
                            latitude >= quantile(latitude, 0.75) ~
                              "North (upper quartile by latitude)",
                            TRUE ~ "Middle") %>%
           factor(levels=c("South (lower quartile by latitude)",
                           # "Middle",
                           "North (upper quartile by latitude)"))) %>%
  na.omit() %>%
  mutate(positivity = positive / total_tested) %>%
  {ggplot(.) +
      aes(x=as.numeric(month), y=total_tested) +
      geom_quasirandom(aes(group=month), varwidth=TRUE, size=0.25) +
      geom_smooth(colour="#FC7A1E", method="gam", formula=y~s(x, bs="cs")) +
      scale_y_log10(breaks=c(0:5, 10*1:5, 100)) +
      ylab("Tests conducted") +
      facet_wrap(~region, nrow=1) +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      ggplot(.) +
      aes(x=as.numeric(month), y=positivity) +
      geom_quasirandom(aes(group=month), varwidth=TRUE, size=0.25) +
      geom_smooth(colour="#FC7A1E", method="gam", formula=y~s(x, bs="cs")) +
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
predictions %>%
  filter(near(latitude, 28, tol=0.5) | near(latitude, 33, tol=0.5) | near(latitude, 37, tol=0.5) | near(latitude, 43, tol=0.5)) %>%
  mutate(latitude = latitude %>% round() %>% str_c("°") %>% factor(ordered=TRUE)) %>%
  {ggplot(.) +
      aes(x=month, colour=latitude) +
      geom_smooth(aes(y=pred_zi), size=1, se=FALSE) +
      scale_y_continuous(name="Detection", limits=0:1) +
      scale_x_continuous(name=NULL, breaks=1:12, labels=month.name) +
      scale_colour_manual(name="Latitude:", values=c("#FC7A1E", "#00A676", "#934D77", "#485696")) +
      theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1),
            axis.ticks.y=element_blank(),
            legend.position=c(0.725, 0.12),
            legend.background=element_rect(fill=alpha("#FFFFFF", 0.8)),
            legend.direction="horizontal")}
ggsave(filename="../figures/detection_by_month_and_latitude.png", width=7, height=4)

predictions %>%
  filter(near(latitude, 28, tol=0.5) | near(latitude, 33, tol=0.5) | near(latitude, 37, tol=0.5) | near(latitude, 43, tol=0.5)) %>%
  mutate(latitude = latitude %>% round() %>% str_c("°") %>% factor(ordered=TRUE)) %>%
  {ggplot(.) +
      aes(x=month, colour=latitude) +
      geom_smooth(aes(y=pred), size=1, se=FALSE) +
      scale_y_continuous(name="Percent Positivity", limits=0:1, breaks=0:4/4,
                         labels=c("0%", "25%", "50%", "75%", "100%")) +
      scale_x_continuous(name=NULL, breaks=1:12, labels=month.name) +
      scale_colour_manual(name="Latitude:", values=c("#FC7A1E", "#00A676", "#934D77", "#485696")) +
      theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1),
            axis.ticks.y=element_blank(),
            legend.position=c(0.725, 0.12),
            # legend.title=element_blank(),
            legend.background=element_rect(fill=alpha("#FFFFFF", 0.8)),
            legend.direction="horizontal")}
ggsave(filename="../figures/prevalence_by_month_and_latitude.png", width=7, height=4)

#### DENSITY LATITUDE INTERACTION PLOT  # # # # # # # # # # # # # # # # # # ####
predictions %>%
  filter(near(density, sort(unique(density))[c(2,4,6,8)])) %>%
  arrange(desc(density)) %>%
  mutate(density = factor(density, levels=sort(unique(density)), labels=c("Low", "Intermediate", "High", "Very High"))) %>%
  {ggplot(.) +
      aes(x=latitude) +
      geom_smooth(aes(y=pred_zi, colour=density), size=1, se=FALSE) +
      geom_rug(size=0.67, alpha=0.2, data=distinct(rug_data, fips, latitude)) +
      scale_y_continuous(name="Detection", limits=0:1) +
      scale_x_continuous(expand=expansion()) +
      scale_colour_manual(values=c("#485696", "#934D77", "#00A676", "#FC7A1E"),
                          name=expression(paste("Population Density (people / ", km^2, ")"))) +
      theme(plot.margin=margin(b=0),
            legend.position=c(0.76, 0.23),
            legend.direction="vertical",
            axis.title.x=element_blank())}
ggsave(filename="../figures/detection_by_latitude_and_density.png", width=7, height=5)

predictions %>%
  filter(near(density, sort(unique(density))[c(2,4,6,8)])) %>%
  arrange(desc(density)) %>%
  mutate(density = factor(density, levels=sort(unique(density)), labels=c("Low", "Intermediate", "High", "Very High"))) %>%
  {ggplot(.) +
      aes(x=latitude) +
      geom_smooth(aes(y=pred, colour=density), size=1, se=FALSE) +
      geom_rug(size=0.67, alpha=0.2, data=distinct(rug_data, fips, latitude)) +
      scale_y_continuous(name="Percent Positivity", limits=0:1, breaks=0:4/4,
                         labels=c("0%", "25%", "50%", "75%", "100%")) +
      scale_x_continuous(expand=expansion()) +
      scale_colour_manual(values=c("#485696", "#934D77", "#00A676", "#FC7A1E"),
                          name=expression(paste("Population Density (people / ", km^2, ")"))) +
      theme(plot.margin=margin(b=0),
            legend.position=c(0.76, 0.23),
            legend.direction="vertical",
            axis.title.x=element_blank())}
ggsave(filename="../figures/prevalence_by_latitude_and_density.png", width=7, height=5)

#### TEMPORAL LAG AND SPATIAL EFFECT  # # # # # # # # # # # # # # # # # # # ####
predictions %>%
  group_by(spatial_lag, temporal_lag) %>%
  summarise(across(c(pred, pred_zi), mean)) %>%
  pivot_longer(c(spatial_lag, temporal_lag), names_to="variable", values_to="value") %>%
  pivot_longer(c(pred, pred_zi), names_to="prediction", values_to="predicted value") %>%
  ggplot() +
  aes(x=value, y=`predicted value`, colour=variable) +
  geom_smooth(size=1) +
  facet_wrap(~prediction, ncol=1, scales="free") +
  scale_colour_manual(values=c("#FC7A1E", "#00A676", "#934D77", "#485696")) +
  theme(axis.ticks.y=element_blank(),
        legend.position=c(0.68, 0.63),
        legend.title=element_blank(),
        legend.background=element_rect(fill=alpha("#FFFFFF", 0.8)),
        legend.direction="horizontal")

#### HABITAT (NOT SIGNIFICANT)  # # # # # # # # # # # # # # # # # # # # # # ####
cutoff_detection <- 0.16; cutoff_positivity <- 0.17; cutoff_persistence <- 0.16;
(monthly_data %>% mutate(detection=as.integer((positive / total_tested) > 0)) %>%
   {ggplot(.) +
       aes(x=habitat, y=detection) +
       geom_quasirandom(groupOnX=FALSE, size=0.25, varwidth=TRUE, width=0.1, alpha=0.05) +
       # geom_smooth(size=1, colour="#FC7A1E") +
       geom_smooth(size=1, colour="#FC7A1E", formula=y~x,
                   method="lm", data=filter(., habitat <= cutoff_detection)) +
       stat_cor(data=filter(., habitat <= cutoff_detection), colour="#FC7A1E", label.x=0.5, label.y=0.4,
                aes(label=paste(after_stat(rr.label), after_stat(p.label), sep="~`,`~")),
                method="spearman") +
       geom_smooth(size=1, colour="#485696", formula=y~x,
                   method="lm", data=filter(., habitat > cutoff_detection)) +
       stat_cor(data=filter(., habitat > cutoff_detection), colour="#485696", label.x=0.5, label.y=0.25,
                aes(label=paste(after_stat(rr.label), after_stat(p.label), sep="~`,`~")),
                method="spearman") +
       scale_y_continuous(name="Detection", breaks=0:1)}) +

  (monthly_data %>%
     mutate(positivity = positive / total_tested) %>%
     {ggplot(.) +
         aes(x=habitat, y=positivity) +
         geom_point(size=0.25, alpha=0.05) +
         # geom_smooth(size=1, colour="#FC7A1E") +
         geom_smooth(size=1, colour="#FC7A1E", formula=y~x,
                     method="lm", data=filter(., habitat <= cutoff_positivity)) +
         stat_cor(data=filter(., habitat <= cutoff_positivity), colour="#FC7A1E", label.x=0.5, label.y=0.4,
                  aes(label=paste(after_stat(rr.label), after_stat(p.label), sep="~`,`~")),
                  method="spearman") +
         geom_smooth(size=1, colour="#485696", formula=y~x,
                     method="lm", data=filter(., habitat > cutoff_positivity)) +
         stat_cor(data=filter(., habitat > cutoff_positivity), colour="#485696", label.x=0.5, label.y=0.25,
                  aes(label=paste(after_stat(rr.label), after_stat(p.label), sep="~`,`~")),
                  method="spearman") +
         scale_y_continuous(name="Percent Positivity", breaks=0:4/4,
                            labels=c("0%", "25%", "50%", "75%", "100%"))}) +

  (persistence_data %>%
     mutate(persistence = as.integer(persistence)) %>%
     {ggplot(.) +
         aes(x=habitat, y=persistence) +
         geom_quasirandom(groupOnX=FALSE, size=0.25, varwidth=TRUE, width=0.1, alpha=0.05) +
         # geom_smooth(size=1, colour="#FC7A1E") +
         geom_smooth(size=1, colour="#FC7A1E", formula=y~x,
                     method="lm", data=filter(., habitat <= cutoff_persistence)) +
         stat_cor(data=filter(., habitat <= cutoff_persistence), colour="#FC7A1E", label.x=0.5, label.y=0.4,
                  aes(label=paste(after_stat(rr.label), after_stat(p.label), sep="~`,`~")),
                  method="spearman") +
         geom_smooth(size=1, colour="#485696", formula=y~x,
                     method="lm", data=filter(., habitat > cutoff_persistence)) +
         stat_cor(data=filter(., habitat > cutoff_persistence), colour="#485696", label.x=0.5, label.y=0.25,
                  aes(label=paste(after_stat(rr.label), after_stat(p.label), sep="~`,`~")),
                  method="spearman") +
         scale_y_continuous(name="Persistence", breaks=0:1)}) +

  plot_layout(ncol=1)
ggsave(filename="../figures/habitat_piecewise_regression.png", width=5, height=7)
