######################
### ALL MODEL FITS ###
######################

rm(list=ls())

### Install packages, their dependencies, and specific versions locally to ensure a reproducible package environment
renv::use(
  "BH@1.78.0-0",
  "Bayesrel@0.7.1",
  "Brobdingnag@1.2-7",
  "DBI@1.1.2",
  "DT@0.20",
  "HDInterval@0.2.2",
  "LaplacesDemon@16.1.6",
  "MASS@7.3-55",
  "Matrix@1.4-0",
  "R6@2.5.1",
  "RColorBrewer@1.1-2",
  "Rcpp@1.0.8",
  "RcppArmadillo@0.10.8.1.0",
  "RcppEigen@0.3.3.9.1",
  "RcppParallel@5.1.5",
  "Rdpack@2.1.3",
  "StanHeaders@2.21.0-7",
  "abind@1.4-5",
  "arrayhelpers@1.1-0",
  "askpass@1.1",
  "assertthat@0.2.1",
  "backports@1.4.1",
  "base64enc@0.1-3",
  "bayesplot@1.8.1",
  "bit64@4.0.5",
  "bit@4.0.4",
  "blob@1.2.2",
  "bridgesampling@1.1-2",
  "brms@2.16.3",
  "broom@0.7.12",
  "bslib@0.3.1",
  "cachem@1.0.6",
  "callr@3.7.0",
  "cellranger@1.1.0",
  "checkmate@2.0.0",
  "cli@3.1.1",
  "clipr@0.7.1",
  "coda@0.19-4",
  "colorspace@2.0-2",
  "colourpicker@1.1.1",
  "commonmark@1.7",
  "cpp11@0.4.2",
  "crayon@1.4.2",
  "crosstalk@1.2.0",
  "curl@4.3.2",
  "data.table@1.14.2",
  "dbplyr@2.1.1",
  "desc@1.4.0",
  "digest@0.6.29",
  "distributional@0.3.0",
  "dplyr@1.0.8",
  "dtplyr@1.2.1",
  "dygraphs@1.1.1.6",
  "ellipsis@0.3.2",
  "evaluate@0.14",
  "fansi@1.0.2",
  "farver@2.1.0",
  "fastmap@1.1.0",
  "fontawesome@0.2.2",
  "forcats@0.5.1",
  "fs@1.5.2",
  "future@1.23.0",
  "gargle@1.2.0",
  "generics@0.1.2",
  "ggdist@3.0.1",
  "ggplot2@3.3.5",
  "ggridges@0.5.3",
  "globals@0.14.0",
  "glue@1.6.1",
  "googledrive@2.0.0",
  "googlesheets4@1.0.0",
  "gridExtra@2.3",
  "gtable@0.3.0",
  "gtools@3.9.2",
  "haven@2.4.3",
  "highr@0.9",
  "hms@1.1.1",
  "htmltools@0.5.2",
  "htmlwidgets@1.5.4",
  "httpuv@1.6.5",
  "httr@1.4.2",
  "ids@1.0.1",
  "igraph@1.2.11",
  "inline@0.3.19",
  "isoband@0.2.5",
  "jquerylib@0.1.4",
  "jsonlite@1.7.3",
  "knitr@1.37",
  "labeling@0.4.2",
  "later@1.3.0",
  "lattice@0.20-45",
  "lavaan@0.6-10",
  "lazyeval@0.2.2",
  "lifecycle@1.0.1",
  "listenv@0.8.0",
  "loo@2.4.1",
  "lubridate@1.8.0",
  "magrittr@2.0.2",
  "markdown@1.1",
  "matrixStats@0.61.0",
  "mgcv@1.8-38",
  "mime@0.12",
  "miniUI@0.1.1.1",
  "mnormt@2.0.2",
  "modelr@0.1.8",
  "munsell@0.5.0",
  "mvtnorm@1.1-3",
  "nleqslv@3.3.2",
  "nlme@3.1-155",
  "numDeriv@2016.8-1.1",
  "openssl@1.4.6",
  "packrat@0.7.0",
  "parallelly@1.30.0",
  "pbivnorm@0.6.0",
  "pillar@1.7.0",
  "pkgbuild@1.3.1",
  "pkgconfig@2.0.3",
  "plyr@1.8.6",
  "posterior@1.2.0",
  "prettyunits@1.1.1",
  "processx@3.5.2",
  "progress@1.2.2",
  "promises@1.2.0.1",
  "ps@1.6.0",
  "purrr@0.3.4",
  "rappdirs@0.3.3",
  "rbibutils@2.2.7",
  "readr@2.1.2",
  "readxl@1.3.1",
  "rematch2@2.1.2",
  "rematch@1.0.1",
  "renv@0.15.2",
  "reprex@2.0.1",
  "reshape2@1.4.4",
  "rlang@1.0.1",
  "rmarkdown@2.11",
  "rprojroot@2.0.2",
  "rsconnect@0.8.25",
  "rstan@2.21.3",
  "rstantools@2.1.1",
  "rstudioapi@0.13",
  "rvest@1.0.2",
  "sass@0.4.0",
  "scales@1.1.1",
  "selectr@0.4-2",
  "shiny@1.7.1",
  "shinyjs@2.1.0",
  "shinystan@2.5.0",
  "shinythemes@1.2.0",
  "sourcetools@0.1.7",
  "stringi@1.7.6",
  "stringr@1.4.0",
  "svUnit@1.0.6",
  "sys@3.4",
  "tensorA@0.36.2",
  "threejs@0.3.3",
  "tibble@3.1.6",
  "tidybayes@3.0.2",
  "tidyr@1.2.0",
  "tidyselect@1.1.1",
  "tidyverse@1.3.1",
  "tinytex@0.36",
  "tmvnsim@1.0-2",
  "tzdb@0.2.0",
  "utf8@1.2.2",
  "uuid@1.0-3",
  "vctrs@0.3.8",
  "viridisLite@0.4.0",
  "vroom@1.5.7",
  "withr@2.4.3",
  "xfun@0.29",
  "xml2@1.3.3",
  "xtable@1.8-4",
  "xts@0.12.1",
  "yaml@2.2.2",
  "zoo@1.8-9"
)

memory.limit(size=560000)

# Load libraries
library(brms) # for implementing Bayesian multilevel analysis
library(tidyverse) # wrangling, plotting, etc.
library(tidybayes) # for preparing and visualizing prior and posterior draws
library(bayesplot) # for prior and posterior predictive plots
library(rstan) # for fitting models directly in Stan
library(Bayesrel) # for Bayesian reliability analysis

# renv::embed() # to initialize and embed lock file in script; generates the renv::use(...) above

# RStan settings
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# to load all fitted models
if (file.exists("pandemic_panel_fits.RData")) base::load(file = "pandemic_panel_fits.RData")

# Read in data
data <- read.csv("main_data.csv") # data for main analysis
d <- data

df_rel <- read.csv("reliability_data.csv")

str(df_rel)

##############################################
### Data preparation: Reliability analysis ###
##############################################

### We only perform the analysis on wave 1 data, since some of the other religiosity variables are unavailable for subsequent waves.

### re-score variables
df_rel$y[df_rel$y==3] <- NA # in wave 1, 3 = "neither agree nor disagree"
df_rel$y[df_rel$y==4] <- 3 # re-scoring
df_rel$y[df_rel$y==5] <- 4 # re-scoring
df_rel$y[df_rel$y==6] <- NA # in wave 1, 6 = "don't know"

hist(df_rel$y) # check

### re-score variables
df_rel$ritual[df_rel$ritual == 5] <- NA # here, 5s are "don't know" -- turn into NAs
table(df_rel$ritual) # check

df_rel$church[df_rel$church == 7] <- NA # here, 7s are "don't know"  -- turn into NAs
table(df_rel$church) # check

df_rel$belief.n <- df_rel$belief # create new column with reversed scores, so that low is less religious
df_rel$belief.n[df_rel$belief.n==4] <- NA # 4s are "don't know" -- turn into NAs
df_rel$belief.n <- abs(df_rel$belief.n - 4) # reverse-scoring: "atheist" = 1, "non-believer" = 2, "believer" = 3 
table(df_rel$belief.n) # check variable, compare with original:
table(df_rel$belief)

df_rel$god.n <- df_rel$god # create new column with reversed scores, so that low is less religious
df_rel$god.n[df_rel$god.n==5] <- NA # 5s are "don't know" -- turn into NAs
df_rel$god.n[df_rel$god.n==3] <- NA # 3s are "don't know what to believe" -- turn into NAs
df_rel$god.n[df_rel$god.n==4] <- 3 # so that the re-structured responses options are in order
df_rel$god.n <- abs(df_rel$god.n - 4) # "no god" = 1, "spiritual force" = 2, "personal god = 3
table(df_rel$god.n) # check variable, compare with original:
table(df_rel$god)

df_rel_rel <- df_rel # final reliability dataset
df_rel_rel$id <- NULL # remove id columns
df_rel_rel$god <- NULL # remove original, now reverse-scored column
df_rel_rel$belief <- NULL  # remove original, now reverse-scored column

str(df_rel_rel)

############################
### Reliability analysis ###
############################

set.seed(2021) # for reproducibility

rel_rel <- strel(data = df_rel_rel, item.dropped = TRUE, missing = "pairwise") # conduct reliability analysis, with drop.item analysis
summary(rel_rel) # summarize

set.seed(NULL) # for reproducibility

#############################################
### Data preparation: Multilevel modeling ###
#############################################

### re-scoring outcome variable
d$y[d$wave == 1][d$y[d$wave == 1]==3] <- NA # (in wave 1, 3 = "neither agree nor disagree")
d$y[d$wave == 1][d$y[d$wave == 1]==4] <- 3 # re-scoring
d$y[d$wave == 1][d$y[d$wave == 1]==5] <- 4 # re-scoring
d$y[d$wave == 1][d$y[d$wave == 1]==6] <- NA # (in wave 1, 6 = "don't know")
hist(d$y[d$wave == 1])

d$y[d$wave == 2][d$y[d$wave == 2]==5] <- NA # (in wave 2, 5 = "prefer not to respond")
d$y[d$wave == 3][d$y[d$wave == 3]==5] <- NA # (in wave 3, 5 = "prefer not to respond")
d$y[d$wave == 2] <- abs(d$y[d$wave == 2] - 5) # rescoring so that low responses are less religious
hist(d$y[d$wave == 2])
d$y[d$wave == 3] <- abs(d$y[d$wave == 3] - 5) # rescoring so that low responses are less religious
hist(d$y[d$wave == 3])

hist(d$y)

### re-scoring covariates
# lump together response 2, 3, and 4, since they capture the same education level (years of formal education)
d$edu[d$edu==3] <- 2
d$edu[d$edu==4] <- 2
d$edu[d$edu==5] <- 3
d$edu[d$edu==6] <- 4
d$edu[d$edu==7] <- 5
d$edu[d$edu==8] <- 6

# Propagate "known" values to missing values in edu and gender:
# when edu[wave==1] = NA, assume that edu is the same in wave 1 as in wave 2 and/or 3
d$edu[d$wave==1] <- ifelse(is.na(d$edu[d$wave==1]), d$edu[d$wave==2], d$edu[d$wave==1]) # if edu is missing in wave 1, take edu from wave 2
d$edu[d$wave==1] <- ifelse(is.na(d$edu[d$wave==1]), d$edu[d$wave==3], d$edu[d$wave==1]) # if edu is missing in wave 1, take edu from wave 3
d$edu <- d$edu[d$wave==1] # set educational level at baseline

# when gender[wave==2] = NA, assume that gender is the same in wave 2 as in wave 1
d$gender[d$wave==1] <- ifelse(is.na(d$gender[d$wave==1]), d$gender[d$wave==2], d$gender[d$wave==1])
# when gender[wave==3] = NA, assume that gender is the same in wave 3 as in wave 1
d$gender[d$wave==3] <- ifelse(is.na(d$gender[d$wave==3]), d$gender[d$wave==1], d$gender[d$wave==3])

# set reference values for covariates
d$edu <- factor(d$edu, c(4,5,6,1,2,3), ordered = TRUE) # set to middle response value
d$health <- factor(d$health, c(3,4,5,1,2), ordered = TRUE) # set to middle response value
d$age.c <- with(d, age[d$wave==1] - mean(age[d$wave==1], na.rm = TRUE)) # age at baseline (wave 1) and center age to sample mean 
d$hinc <- factor(d$hinc, c(6,7,8,9,10,11,1,2,3,4,5), ordered = TRUE) # set to middle response value
# gender = 1 = woman

d_full <- data.frame(id = as.factor(d$id), 
                     wave = as.integer(d$wave),
                     y = factor(as.factor(d$y), ordered = TRUE), 
                     health = d$health,
                     gender = as.factor(d$gender),
                     age = as.numeric(d$age),
                     age.c = d$age.c, 
                     edu = d$edu, 
                     hinc = d$hinc)

str(d_full) # check structure of final data set

### Model 0: Unconditional (marginal) change model 

## Prior predictive simulations

# get list of prior parameters
get_prior(y ~ 1 + mo(wave) + (1 + mo(wave) | id),
          data = d_full, 
          family=cumulative("logit"))

# set priors
priors_m0 <- set_prior("normal(0,0.5)", class = "b") +
  set_prior("dirichlet(2)", class="simo", coef = "mowave1") +
  set_prior("normal(0,10)", class = "Intercept") +
  set_prior("exponential(1)", class = "sd") +
  set_prior("lkj_corr_cholesky(4)", class = "cor")

priorcheck_m0 <- brm(y ~ 1 + mo(wave) + (1 + mo(wave) | id),
                     data = d_full,
                     prior = priors_m0,
                     family = cumulative("logit"),
                     chains = 2, iter = 1000,
                     sample_prior = "only",
                     seed = 2021)

summary(priorcheck_m0)

set.seed(2021)

yrep_pc_m0 <- posterior_predict(priorcheck_m0)

priorplot_m0 <- ppc_bars_grouped(y = as.numeric(priorcheck_m0$data[["y"]]),
                yrep = yrep_pc_m0,
                group = priorcheck_m0$data[["wave"]],
                freq = FALSE) # easier to compare groups with proportions, instead of counts

set.seed(NULL)

## Fit Model 0 
m0 <- brm(y ~ 1 + mo(wave) + (1 + mo(wave) | id), # marginal model (listwise deletion)
                  data = d_full,
                  family=cumulative("logit"),
                  prior = priors_m0,
                  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99),
                  seed = 2021)

summary(m0)

set.seed(2021)

yrep_m0 <- posterior_predict(m0)

postplot_m0 <- ppc_bars_grouped(y = as.numeric(m0$data[["y"]]),
                  yrep = yrep_m0,
                  group = m0$data[["wave"]],
                  freq = FALSE) # easier to compare groups with proportions, instead of counts

set.seed(NULL)

## Fit Model 0_complete
d_complete <- d_full[complete.cases(d_full$y),] %>% group_by(id) %>% filter(n()>2) %>% as.data.frame() # include only participants who completed all three waves

m0_complete <- brm(y ~ 1 + mo(wave) + (1 + mo(wave) | id),
                             data = d_complete,
                             family=cumulative("logit"),
                             prior = priors_m0,
                             cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99),
                             seed = 2021)

summary(m0_complete)

set.seed(2021)

yrep_m0_complete <- posterior_predict(m0_complete)

postplot_m0_complete <- ppc_bars_grouped(y = as.numeric(m0_complete$data[["y"]]),
                  yrep = yrep_m0_complete,
                  group = m0_complete$data[["wave"]],
                  freq = FALSE)

set.seed(NULL)

### Model 1

## Prior predictive simulations

# get list of prior parameters
get_prior(y ~ 1 + mo(wave) + mo(edu) + mo(hinc) + mo(health) + age.c + gender + mo(wave)*mo(hinc) + mo(wave)*mo(health) + (1 + mo(wave) | id),
          data = d_full, 
          family=cumulative("logit"))

# set priors
priors_m1 <- set_prior("normal(0,0.5)", class = "b") +
  set_prior("normal(0,10)", class = "Intercept") +
  set_prior("exponential(1)", class = "sd") +
  set_prior("lkj_corr_cholesky(4)", class = "cor") +
  set_prior("dirichlet(2)", class="simo", 
            coef = c("mowave1","moedu1","mohealth1",
                     "mowave:mohealth1","mowave:mohealth2",
                     "mohinc1", "mowave:mohinc1","mowave:mohinc2"))
    
priorcheck_m1 <- brm(y ~ 1 + mo(wave) + mo(edu) + mo(hinc) + mo(health) + age.c + gender + mo(wave)*mo(hinc) + mo(wave)*mo(health) + (1 + mo(wave) | id),
                     data = d_full,
                     prior = priors_m1,
                     family = cumulative("logit"),
                     chains = 2, iter = 1000,
                     sample_prior = "only",
                     seed = 2021)

summary(priorcheck_m1)

set.seed(2021)

yrep_pc_m1 <- posterior_predict(priorcheck_m1)

priorplot_m1 <- ppc_bars_grouped(y = as.numeric(priorcheck_m1$data[["y"]]),
                  yrep = yrep_pc_m1,
                  group = priorcheck_m1$data[["wave"]],
                  freq = FALSE) # easier to compare groups with proportions, instead of counts

set.seed(NULL)

### Fit Model 1: Change model including full covariate set (listwise deletion)
m1 <- brm(y ~ 1 + mo(wave) + mo(edu) + mo(hinc) + mo(health) + age.c + gender + mo(wave)*mo(hinc) + mo(wave)*mo(health) + (1 + mo(wave) | id),
                  data = d_full, 
                  family=cumulative("logit"),
                  prior = priors_m1,
                  cores = 4, chains = 4, iter = 6000, control = list(adapt_delta = 0.99),
                  seed = 2021)

summary(m1)

set.seed(2021)

yrep_m1 <- posterior_predict(m1)

postplot_m1 <- ppc_bars_grouped(y = as.numeric(m1$data[["y"]]),
                  yrep = yrep_m1,
                  group = m1$data[["wave"]],
                  freq = FALSE)

set.seed(NULL)

### Fit Model 1_complete
d_full_complete <- d_full[complete.cases(d_full),] %>% group_by(id) %>% filter(n()>2) %>% ungroup() %>% as.data.frame() # include only participants who completed all three waves

m1_complete <- brm(y ~ 1 + mo(wave) + mo(edu) + mo(hinc) + mo(health) + age.c + gender + mo(wave)*mo(hinc) + mo(wave)*mo(health) + (1 + mo(wave) | id),
                   data = d_full_complete, 
                   family=cumulative("logit"),
                   prior = priors_m1,
                   cores = 4, chains = 4, iter = 6000, control = list(adapt_delta = 0.99),
                   seed = 2021)

summary(m1_complete)

set.seed(2021)

yrep_m1_complete <- posterior_predict(m1_complete)

postplot_m1_complete <- ppc_bars_grouped(y = as.numeric(m1_complete$data[["y"]]),
                  yrep = yrep_m1_complete,
                  group = m1_complete$data[["wave"]],
                  freq = FALSE)

set.seed(NULL)

### "Unconditional" interaction models

## health model: prior predictive check

# set priors
priors_health <- set_prior("normal(0,0.5)", class = "b") +
  set_prior("normal(0,10)", class = "Intercept") +
  set_prior("exponential(1)", class = "sd") +
  set_prior("lkj_corr_cholesky(4)", class = "cor") +
  set_prior("dirichlet(2)", class="simo", coef = "mohealth1") +
  set_prior("dirichlet(2)", class="simo", coef = "mowave1") +
  set_prior("dirichlet(2)", class="simo", coef = "mowave:mohealth1") + 
  set_prior("dirichlet(2)", class="simo", coef = "mowave:mohealth2")

# sample from the priors
pc_health_model <- brm(y ~ 1 + mo(wave) + mo(wave)*mo(health) + (1 + mo(wave) + mo(wave)*mo(health) | id), 
                       data = d_full,
                       family = cumulative("logit"), 
                       prior = priors_health,
                       chains = 2, iter = 1000,
                       seed = 2021,
                       sample_prior = "only")

summary(pc_health_model)

set.seed(2021)

yrep_pc_health <- posterior_predict(pc_health_model)

priorplot_health <- ppc_bars_grouped(y = as.numeric(pc_health_model$data[["y"]]),
                  yrep = yrep_pc_health,
                  group = pc_health_model$data[["wave"]],
                  freq = FALSE)

set.seed(NULL)

### "Unconditional" interaction models

## income model: prior predictive check

# set priors
priors_hinc <- set_prior("normal(0,0.5)", class = "b") +
  set_prior("normal(0,10)", class = "Intercept") +
  set_prior("exponential(1)", class = "sd") +
  set_prior("lkj_corr_cholesky(4)", class = "cor") +
  set_prior("dirichlet(2)", class="simo", coef = "mowave1") +
  set_prior("dirichlet(2)", class="simo", coef = "mohinc1") +
  set_prior("dirichlet(2)", class="simo", coef = "mowave:mohinc1") +
  set_prior("dirichlet(2)", class="simo", coef = "mowave:mohinc2")

# sample from the priors
pc_hinc_model <- brm(y ~ 1 + mo(wave) + mo(wave)*mo(hinc) + (1 + mo(wave) + mo(wave)*mo(hinc) | id), 
                       data = d_full,
                       family = cumulative("logit"), 
                       prior = priors_hinc,
                       chains = 2, iter = 1000,
                       seed = 2021,
                       sample_prior = "only")

summary(pc_hinc_model)

set.seed(2021)

yrep_pc_hinc <- posterior_predict(pc_hinc_model)

priorplot_hinc <- ppc_bars_grouped(y = as.numeric(pc_hinc_model$data[["y"]]),
                  yrep = yrep_pc_hinc,
                  group = pc_hinc_model$data[["wave"]],
                  freq = FALSE)

set.seed(NULL)

# fit model
health_model <- brm(y ~ 1 + mo(wave) + mo(wave)*mo(health) + (1 + mo(wave) + mo(wave)*mo(health) | id), 
                        data = d_full,
                        family = cumulative("logit"), 
                        prior = priors_health,
                        seed = 2021,
                        cores = 4, chains = 4, iter = 6000, control = list(adapt_delta = 0.99))

summary(health_model)

set.seed(2021)

yrep_health <- posterior_predict(pc_health_model)

postplot_health <- ppc_bars_grouped(y = as.numeric(pc_health_model$data[["y"]]),
                  yrep = yrep_health,
                  group = pc_health_model$data[["wave"]],
                  freq = FALSE)

set.seed(NULL)

# fit model
hinc_model <- brm(y ~ 1 + mo(wave) + mo(wave)*mo(hinc) + (1 + mo(wave) + mo(wave)*mo(hinc) | id), 
                        data = d_full,
                        family = cumulative("logit"), 
                        prior = priors_hinc,
                        seed = 2021,
                        cores = 4, chains = 4, iter = 6000, control = list(adapt_delta = 0.99))

summary(hinc_model)

set.seed(2021)

yrep_hinc <- posterior_predict(hinc_model)

postplot_hinc <- ppc_bars_grouped(y = as.numeric(hinc_model$data[["y"]]),
                  yrep = yrep_hinc,
                  group = hinc_model$data[["wave"]],
                  freq = FALSE)

set.seed(NULL)

### Poststratification (MRP) model

## Select variables to include in MRP model
data_mrp <- with(d_full, 
                 data.frame(id = id, 
                            wave = wave, 
                            y = y, 
                            gender = gender, 
                            age = age[wave==1], # set age to baseline age
                            edu = as.numeric(edu))) # already set to baseline educational level

## Bin covariates according to external data

# to achieve maximal coverage of age groups while retaining reasonable (i.e., not too large, not too small) bin sizes
# consistent with externally available age ranges, we chose the following ranges:

data_mrp$age[data_mrp$age <= 24] <- 1 # below 24 (external age bin: 15-24)
data_mrp$age[data_mrp$age >= 25 & data_mrp$age <= 49] <- 2 # 25-49
data_mrp$age[data_mrp$age >= 50 & data_mrp$age <= 74] <- 3 # 50-74
data_mrp$age[data_mrp$age > 74] <- 4 # above 74; will not be used to poststratify due to lack of external data for this group

data_mrp$age <- as.factor(data_mrp$age) # convert to factor

data_mrp$edu[data_mrp$edu == 1] <- 1 # primary to lower secondary
data_mrp$edu[data_mrp$edu == 2 & data_mrp$edu ==3] <- 2 # upper secondary and post-secondary non-tertiary education
data_mrp$edu[data_mrp$edu >= 4] <- 3 # tertiary education

data_mrp$edu <- as.factor(data_mrp$edu) # convert to factor

d_mrp <- data_mrp[complete.cases(data_mrp$y),] # listwise deletion of missing outcome

str(d_mrp)

## MRP Model

# Set priors
priors_mrp <- set_prior("normal(0,0.5)", class = "b") +
  set_prior("normal(0,10)", class = "Intercept") +
  set_prior("exponential(1)", class = "sd") +
  set_prior("dirichlet(2)", class="simo", coef = "mowave1") +
  set_prior("lkj_corr_cholesky(4)", class = "cor")

# MRP model prior check
priorcheck_m_mrp <- brm(y ~ 1 + mo(wave) + (1 | edu) + (1 | age) + (1 | gender) + (1 + mo(wave) | id),
                 data = d_mrp, 
                 family = cumulative("logit"),
                 prior = priors_mrp,
                 sample_prior = "only",
                 cores = 4, chains = 4, iter = 1000,
                 seed = 2021)

# MRP model with listwise deletion
m_mrp <- brm(y ~ 1 + mo(wave) + (1 | edu) + (1 | age) + (1 | gender) + (1 + mo(wave) | id),
                       data = d_mrp, 
                       family = cumulative("logit"),
                       prior = priors_mrp,
                       cores = 4, chains = 4, iter = 8000, control = list(adapt_delta = 0.99),
                       seed = 2021)

set.seed(2021)

yrep_mrp <- posterior_predict(m_mrp)

postplot_mrp <- ppc_bars_grouped(y = as.numeric(m_mrp$data[["y"]]),
                  yrep = yrep_mrp,
                  group = m_mrp$data[["wave"]],
                  freq = FALSE)

set.seed(NULL)

# poststratification takes places in script for supplementary plots

## Imputation model

# create a dataset with complete outcome; treat missing categorical covariates as numeric to allow imputation.....
y_notmiss <- which( !is.na(d_full$y)) 
d_imp <- with(d_full, 
              data.frame(
              id = id[y_notmiss],
              wave = wave[y_notmiss],
              y = y[y_notmiss],
              edu = as.numeric(edu[y_notmiss]),
              hinc = as.numeric(hinc[y_notmiss]),
              health = as.numeric(health[y_notmiss]),
              age.c = age.c[y_notmiss],
              gender = gender[y_notmiss]))

str(d_imp)

# specify imputation model formula
imp_formula <- bf(y ~ 1 + mo(wave) + mi(edu) + mi(hinc) + mi(health) + age.c + gender + mo(wave)*mi(hinc) + mo(wave)*mi(health) + (1 + mo(wave) | id), family=cumulative("logit")) +
               bf(edu | mi() ~ age.c + gender, family=gaussian) + 
               bf(hinc | mi() ~ mo(wave) + age.c + gender + mi(edu) + mi(health), family=gaussian) +
               bf(health | mi() ~ mo(wave) + age.c + gender + mi(edu) + mi(hinc), family=gaussian) + set_rescor(FALSE)

# Set priors
get_prior(formula = imp_formula, 
          data = d_imp)

priors_imp <- set_prior("normal(0,0.5)", class = "b", resp = c("y", "edu", "health", "hinc")) +
              set_prior("normal(0,10)", class = "Intercept", resp = "y") +
              set_prior("exponential(1)", class = "sd", resp = "y") +
              set_prior("lkj_corr_cholesky(4)", class = "cor") +
              set_prior("dirichlet(2)", class = "simo", coef = c("mowave1", "mowave:mihealth1", "mowave:mihinc1"), resp = "y") +
              set_prior("dirichlet(2)", class = "simo", coef = "mowave1", resp = "health") + 
              set_prior("dirichlet(2)", class = "simo", coef = "mowave1", resp = "hinc")

# prior check
m1_imp_ppc <- brm(formula = imp_formula, 
              data = d_imp,
              prior = priors_imp,
              sample = "only",
              iter = 1000, chains = 4, cores = 4,
              seed = 2021)

prior_summary(m1_imp_ppc)

set.seed(2021)

yrep_pc_imp <- posterior_predict(m1_imp_ppc, resp = "y")

priorplot_imp <- ppc_bars_grouped(y = as.numeric(m1_imp_ppc$data[["y"]]),
                  yrep = yrep_pc_imp,
                  group = m1_imp_ppc$data[["wave"]],
                  resp = "y",
                  freq = FALSE)

set.seed(NULL)

# we want to constrain the imputed values to realistic ranges; however, this cannot (AFAIK) be specified directly in brms,
# so we have to workaround that by 1) creating an "empty" brmsfit, 2) amend the raw Stan code with the constraints, 
# 3) then feed the amended Stan code (including the data) back into our brmsfit via rstan,
# as suggested by B<fc>rkner here: https://discourse.mc-stan.org/t/creating-a-brmsfit-object-with-a-modified-brms-generated-stan-model/6840

# create "empty" (0 chains) brmsfit object
m1_imp <- brm(formula = imp_formula, 
                 data = d_imp,
                prior = priors_imp,
               chains = 0)

# extract raw Stan code
imp_stancode <- stancode(m1_imp)

# amend Stan raw code: set constraints for imputed values
imp_stancode_edu <- gsub("vector[Nmi_edu] Ymi_edu;", "vector<lower=1,upper=6>[Nmi_edu] Ymi_edu;", imp_stancode, fixed=TRUE)
imp_stancode_hinc <- gsub("vector[Nmi_hinc] Ymi_hinc;", "vector<lower=1,upper=11>[Nmi_hinc] Ymi_hinc;", imp_stancode_edu, fixed=TRUE)
imp_stancode_health <- gsub("vector[Nmi_health] Ymi_health;", "vector<lower=1,upper=5>[Nmi_health] Ymi_health;", imp_stancode_hinc, fixed=TRUE)
imp_stancode <- imp_stancode_health

# create standata
imp_standata <- make_standata(imp_formula, d_imp, priors = priors_imp)

# create "empty" (0 chains) stanfit object
imp_stanfit <- rstan::stan(model_code = imp_stancode, data = imp_standata, chains = 0)

# feed the stanfit to the brmsfit
m1_imp$fit <- imp_stanfit

# fit full imputation model via Stan to "empty" brmsfit object
m1_imp <- stats::update(m1_imp, recompile = FALSE, iter = 6000, chains = 4, seed = 2021)

set.seed(2021)

yrep_imp <- posterior_predict(m1_imp, resp = "y")

postplot_imp <- ppc_bars_grouped(y = as.numeric(m1_imp$data[["y"]]),
                  yrep = yrep_imp,
                  group = m1_imp$data[["wave"]],
                  resp = "y",
                  freq = FALSE)

set.seed(NULL)

### Exposure models

# read in exposure data
df_exp <- read.csv("exposure_data.csv")

### re-scoring outcome variable
df_exp$y[df_exp$wave == 1][df_exp$y[df_exp$wave == 1]==3] <- NA # (in wave 1, 3 = "neither agree nor disagree")
df_exp$y[df_exp$wave == 1][df_exp$y[df_exp$wave == 1]==4] <- 3 # re-scoring
df_exp$y[df_exp$wave == 1][df_exp$y[df_exp$wave == 1]==5] <- 4 # re-scoring
df_exp$y[df_exp$wave == 1][df_exp$y[df_exp$wave == 1]==6] <- NA # (in wave 1, 6 = "don't know")

df_exp$y[df_exp$wave == 3][df_exp$y[df_exp$wave == 3]==5] <- NA # (in wave 3, 5 = "prefer not to respond")
df_exp$y[df_exp$wave == 3] <- abs(df_exp$y[df_exp$wave == 3] - 5) # rescoring so that low responses are less religious

df_exp$y <- factor(as.factor(df_exp$y), ordered = TRUE)

str(df_exp)

# Set priors
get_prior(y ~ 1 + mo(wave)*mo(exp_self) + (1 + mo(wave)*mo(exp_self) | id), 
          data = df_exp)

priors_exp <- set_prior("normal(0,0.5)", class = "b") +
              set_prior("normal(0,10)", class = "Intercept") +
              set_prior("exponential(1)", class = "sd") +
              set_prior("lkj_corr_cholesky(4)", class = "cor") +
              set_prior("dirichlet(2)", class="simo", coef = "mowave1")

# "Self model"
priors_exp_self <- priors_exp  +
              set_prior("dirichlet(2)", class="simo", coef = "mowave:moexp_self1") +
              set_prior("dirichlet(2)", class="simo", coef = "mowave:moexp_self2") +
              set_prior("dirichlet(2)", class="simo", coef = "moexp_self1")

m_self <- brm(y ~ 1 + mo(wave)*mo(exp_self) + (1 + mo(wave)*mo(exp_self) | id),
                    data = df_exp,
                    family=cumulative("logit"),
                    prior = priors_exp_self,
                    cores = 4, chains = 4, iter = 6000, control = list(adapt_delta = 0.99),
                    seed = 2021)

# "Household model"
priors_exp_household <- priors_exp  +
  set_prior("dirichlet(2)", class="simo", coef = "mowave:moexp_household1") +
  set_prior("dirichlet(2)", class="simo", coef = "mowave:moexp_household2") +
  set_prior("dirichlet(2)", class="simo", coef = "moexp_household1")

m_household <- brm(y ~ 1 + mo(wave)*mo(exp_household) + (1 + mo(wave)*mo(exp_household) | id),
                    data = df_exp,
                    family=cumulative("logit"),
                    prior = priors_exp_household,
                    cores = 4, chains = 4, iter = 6000, control = list(adapt_delta = 0.99),
                    seed = 2021)

# "Family model"
priors_exp_family <- priors_exp  +
  set_prior("dirichlet(2)", class="simo", coef = "mowave:moexp_family1") +
  set_prior("dirichlet(2)", class="simo", coef = "mowave:moexp_family2") +
  set_prior("dirichlet(2)", class="simo", coef = "moexp_family1")

m_family <- brm(y ~ 1 + mo(wave)*mo(exp_family) + (1 + mo(wave)*mo(exp_family) | id),
                    data = df_exp,
                    family=cumulative("logit"),
                    prior = priors_exp_family,    
                    cores = 4, chains = 4, iter = 6000, control = list(adapt_delta = 0.99),
                    seed = 2021)

# "Relations model"
priors_exp_relations <- priors_exp  +
  set_prior("dirichlet(2)", class="simo", coef = "mowave:moexp_relations1") +
  set_prior("dirichlet(2)", class="simo", coef = "mowave:moexp_relations2") +
  set_prior("dirichlet(2)", class="simo", coef = "moexp_relations1")

m_relations<- brm(y ~ 1 + mo(wave)*mo(exp_relations) + (1 + mo(wave)*mo(exp_relations) | id),
                    data = df_exp,
                    family=cumulative("logit"),
                    prior = priors_exp_relations,  
                    cores = 4, chains = 4, iter = 6000, control = list(adapt_delta = 0.99),
                    seed = 2021)

save.image("pandemic_panel_fits.RData")

### END
