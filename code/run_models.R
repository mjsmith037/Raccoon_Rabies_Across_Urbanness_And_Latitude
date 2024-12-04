library(magrittr)
library(lme4)
library(GLMMadaptive)
library(splines)
library(rlang)
library(tidyverse)

################################################################################
# This script contains code to re-run all three models reported in the
# publication. If you want to run a subset, just comment/omit running the
# code blocks starting with `assign` and `save` for the models you don't need
################################################################################

# get the data (hypothetical data provided for demonstrative purposes)
# NOTE: use of the hypothetical data will result in non-convergence of the models
monthly_data <- read_csv("../data/hypothetical_data.csv", col_types="fdddddddd") %>%
  mutate(latitude = scale(latitude),
         density = scale(density))

# best fitting spline parameters (differ by model type)
df_density <- 1
degree_density <- 1

df_month <- 6
degree_month <- 1

df_latitude <- 5
degree_latitude <- 1

# construct the base model formula used for the positivity and presence models
base_formula <- formula(str_glue("~ bs(month, df={df_month}, degree={degree_month}) +
  bs(latitude, df={df_latitude}, degree={degree_latitude}) +
  bs(density, df={df_density}, degree={degree_density}) +
  bs(density, df={df_density}, degree={degree_density}) : bs(latitude, df={df_latitude}, degree={degree_latitude}) +
  bs(month, df={df_month}, degree={degree_month}) : bs(latitude, df={df_latitude}, degree={degree_latitude}) +
  habitat + temporal_lag + spatial_lag"))

# this will be used for object and file naming purposes. Note we name the splines
# based on the number of parts (df - degree + 1) rather than the df and degree explicitly
model_base <- str_glue("density_{df_density-degree_density+1}_{degree_density}_",
                       "month_{df_month-degree_month+1}_{degree_month}_",
                       "latitude_{df_latitude-degree_latitude+1}_{degree_latitude}")

#### Rabies Percent Positivity Model ####
assign(str_glue("positivity_{model_base}"),
       mixed_model(
         fixed     = update.formula(base_formula, positivity ~ . + log(total_tested)),
         random    = ~ 1 | fips,
         family    = beta.fam(),
         max_coef_value = 1000,
         data      = monthly_data %>%
           mutate(positivity = positive %>%
                    divide_by(total_tested) %>%
                    # we do a slight transformation to adjust
                    # edge cases, avoiding an error when fitting
                    multiply_by(n() - 1) %>% add(0.5) %>% divide_by(n()))
       ))
# save(list=c(str_glue("positivity_{model_base}"), "base_formula"),
#      file=str_glue("../results/positivity_{model_base}.RData"))

#### Rabies Presence Model ####
assign(str_glue("presence_{model_base}"),
       mixed_model(
         fixed     = update.formula(base_formula, presence ~ . + log(total_tested)),
         random    = ~ 1 | fips,
         family    = binomial(),
         data      = monthly_data %>% mutate(presence = (positive > 0) %>% as.factor())
       ))
# save(list=c(str_glue("presence_{model_base}"), "base_formula"),
#      file=str_glue("../results/presence_{model_base}.RData"))

#### Rabies Persistence Model ####
# This model requires restructured data. First define the lag and cutoff parameters
num_months_with_pos <- 5
lag_duration_months <- 11

# transform the original data to identify cases of rabies persistence
# this is not very fast/efficient, so we load if available and save the output
# as a new csv file if we do have to calculate it
if (file.exists(str_glue("../data/persistence_{num_months_with_pos}.{lag_duration_months}_data.csv"))) {
  persistence_data <- read_csv(str_glue("../data/persistence_{num_months_with_pos}.{lag_duration_months}_data.csv"))
} else {
  # lag function adapted from:
  # https://purrple.cat/blog/2018/03/02/multiple-lags-with-tidy-evaluation/
  lags <- function(var, lag_duration){
    var <- enquo(var)
    indices <- seq_len(lag_duration)
    map(indices, ~quo(lag(!!var, !!.x))) %>%
      set_names(sprintf("lag_%s_%02d", quo_text(var), indices))
  }
  persistence_data <- monthly_data %>%
    mutate(year = year - min(year),
           across(c(density, latitude, spatial_lag, temporal_lag), scale)) %>%
    group_by(fips) %>%
    # get the lagged positives and total tested
    mutate(!!!lags(positive, lag_duration_months), !!!lags(total_tested, lag_duration_months)) %>%
    # which lag intervals have positives?
    mutate(across(starts_with("lag_positive"), ~is_greater_than(., 0))) %>%
    rowwise() %>%
    # identify which counties meet the criteria (and tally the total tests over the lag)
    mutate(persistence = sum(c_across(starts_with("lag_positive"))) %>%
             is_weakly_greater_than(num_months_with_pos),
           total_tested_over_lag = sum(c_across(starts_with("lag_total_tested"))),
           .keep="unused") %>%
    # clean up
    ungroup() %>% na.omit() %>% mutate(across(where(is.matrix), ~.[,1])) %>%
    write_csv(str_glue("../data/persistence_{num_months_with_pos}.{lag_duration_months}_data.csv"))
}

assign(str_glue("persistence_{num_months_with_pos}.{lag_duration_months}_{model_base}"),
       mixed_model(
         fixed          = base_formula,
         random         = ~ 1 | fips,
         family         = binomial(),
         max_coef_value = 1000,
         iter_EM        = 0,
         data           = persistence_data
       ))
# save(list=c(str_glue("persistence_{num_months_with_pos}.{lag_duration_months}_{model_base}"), "base_formula"),
#      file=str_glue("../results/persistence_{num_months_with_pos}.{lag_duration_months}_{model_base}.RData"))


