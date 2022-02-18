###########################
### SUPPLEMENTARY PLOTS ###
###########################

### Install packages, their dependencies, and specific versions locally to ensure a reproducible package environment
renv::use(
  "BH@1.75.0-0",
  "Brobdingnag@1.2-6",
  "DBI@1.1.2",
  "DT@0.20",
  "HDInterval@0.2.2",
  "MASS@7.3-54",
  "Matrix@1.3-4",
  "R6@2.5.1",
  "RColorBrewer@1.1-2",
  "Rcpp@1.0.8",
  "RcppEigen@0.3.3.9.1",
  "RcppParallel@5.1.5",
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
  "bslib@0.3.1",
  "cachem@1.0.6",
  "callr@3.7.0",
  "cellranger@1.1.0",
  "checkmate@2.0.0",
  "cli@3.1.1",
  "clipr@0.7.1",
  "coda@0.19-4",
  "codetools@0.2-18",
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
  "dplyr@1.0.7",
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
  "modelr@0.1.8",
  "munsell@0.5.0",
  "mvtnorm@1.1-3",
  "nleqslv@3.3.2",
  "nlme@3.1-153",
  "numDeriv@2016.8-1.1",
  "openssl@1.4.5",
  "packrat@0.7.0",
  "parallelly@1.29.0",
  "patchwork@1.1.1",
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
  "readr@2.1.1",
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
  "tidymodels/broom@HEAD",
  "tidyr@1.2.0",
  "tidyselect@1.1.1",
  "tidyverse@1.3.1",
  "tinytex@0.36",
  "tzdb@0.2.0",
  "utf8@1.2.2",
  "uuid@1.0-3",
  "vctrs@0.3.8",
  "viridis@0.6.2",
  "viridisLite@0.4.0",
  "vroom@1.5.7",
  "withr@2.4.3",
  "xfun@0.29",
  "xml2@1.3.2",
  "xtable@1.8-4",
  "xts@0.12.1",
  "yaml@2.2.2",
  "zoo@1.8-9"
)

# renv::embed() # to initialize and embed lock file in script; generates the renv::use(...) above

### Load packages
library(brms) # for implementing Bayesian multilevel analysis
library(tidyverse) # wrangling, plotting, etc.
library(tidybayes) # for preparing and visualizing prior and posterior draws
library(bayesplot) # for prior and posterior predictive plots
library(patchwork) # for panel plots
library(modelr) # for plotting
library(xtable) # for printing TeX tables

# load all fitted models
if (file.exists("pandemic_panel_fits.RData")) base::load(file = "pandemic_panel_fits.RData")

memory.limit(size=560000)

### PRIOR VS. POSTERIOR PREDICTIVE CHECKS

## m0:

# prior predictive check
set.seed(2021)

yrep_pc_m0 <- posterior_predict(priorcheck_m0)

priorplot_m0 <- ppc_bars_grouped(y = as.numeric(priorcheck_m0$data[["y"]]),
                                 yrep = yrep_pc_m0,
                                 group = priorcheck_m0$data[["wave"]],
                                 prob = 0.95,
                                 size = 0.5,
                                 fatten = 1,
                                 freq = FALSE) + # easier to compare groups (i.e., wave) with proportions, instead of counts
                                 ylim(c(0,1)) + theme(legend.position = "null") + ggtitle("Prior predictive check") +
                                 scale_x_continuous(name="Response options")
set.seed(NULL)

# posterior predictive check
set.seed(2021)

yrep_m0 <- posterior_predict(m0)

postplot_m0 <- ppc_bars_grouped(y = as.numeric(m0$data[["y"]]),
                                yrep = yrep_m0,
                                group = m0$data[["wave"]],
                                prob = 0.95,
                                size = 0.5,
                                fatten = 1,
                                freq = FALSE) + # easier to compare groups (i.e., wave) with proportions, instead of counts
                                ylim(c(0,1)) + theme(legend.position = "null") + ggtitle("Posterior predictive check") +
                                scale_x_continuous(name="Response options")
set.seed(NULL)

# plot
cairo_pdf("m0_ppc.pdf", width = 7, height = 5.5) # start print to pdf

(priorplot_m0 + postplot_m0 + patchwork::plot_layout(ncol = 1, nrow = 2))

dev.off() # end print to pdf

## m1:

# prior predictive check
set.seed(2021)

yrep_pc_m1 <- posterior_predict(priorcheck_m1)

priorplot_m1 <- ppc_bars_grouped(y = as.numeric(priorcheck_m1$data[["y"]]),
                                 yrep = yrep_pc_m1,
                                 group = priorcheck_m1$data[["wave"]],
                                 prob = 0.95,
                                 size = 0.5,
                                 fatten = 1,
                                 freq = FALSE) + # easier to compare groups (i.e., wave) with proportions, instead of counts
  ylim(c(0,1)) + theme(legend.position = "null") + ggtitle("Prior predictive check") +
  scale_x_continuous(name="Response options")
set.seed(NULL)

# posterior predictive check
set.seed(2021)

yrep_m1 <- posterior_predict(m1)

postplot_m1 <- ppc_bars_grouped(y = as.numeric(m1$data[["y"]]),
                                yrep = yrep_m1,
                                group = m1$data[["wave"]],
                                prob = 0.95,
                                size = 0.5,
                                fatten = 1,
                                freq = FALSE) + # easier to compare groups (i.e., wave) with proportions, instead of counts
  ylim(c(0,1)) + theme(legend.position = "null") + ggtitle("Posterior predictive check") +
  scale_x_continuous(name="Response options")
set.seed(NULL)

# plot
cairo_pdf("m1_ppc.pdf", width = 7, height = 5.5) # start print to pdf

(priorplot_m1 + postplot_m1 + patchwork::plot_layout(ncol = 1, nrow = 2))

dev.off() # end print to pdf

## hinc_model:

# prior predictive check
set.seed(2021)

yrep_pc_hinc_model <- posterior_predict(pc_hinc_model)

priorplot_hinc_model <- ppc_bars_grouped(y = as.numeric(pc_hinc_model$data[["y"]]),
                                 yrep = yrep_pc_hinc_model,
                                 group = pc_hinc_model$data[["wave"]],
                                 prob = 0.95,
                                 size = 0.5,
                                 fatten = 1,
                                 freq = FALSE) + # easier to compare groups (i.e., wave) with proportions, instead of counts
  ylim(c(0,1)) + theme(legend.position = "null") + ggtitle("Prior predictive check") +
  scale_x_continuous(name="Response options")
set.seed(NULL)

# posterior predictive check
set.seed(2021)

yrep_hinc_model <- posterior_predict(hinc_model)

postplot_hinc_model <- ppc_bars_grouped(y = as.numeric(hinc_model$data[["y"]]),
                                yrep = yrep_hinc_model,
                                group = hinc_model$data[["wave"]],
                                prob = 0.95,
                                size = 0.5,
                                fatten = 1,
                                freq = FALSE) + # easier to compare groups (i.e., wave) with proportions, instead of counts
  ylim(c(0,1)) + theme(legend.position = "null") + ggtitle("Posterior predictive check") +
  scale_x_continuous(name="Response options")
set.seed(NULL)

# plot
cairo_pdf("hinc_model_ppc.pdf", width = 7, height = 5.5) # start print to pdf

(priorplot_hinc_model + postplot_hinc_model + patchwork::plot_layout(ncol = 1, nrow = 2))

dev.off() # end print to pdf

## health_model:

# prior predictive check
set.seed(2021)

yrep_pc_health_model <- posterior_predict(pc_health_model)

priorplot_health_model <- ppc_bars_grouped(y = as.numeric(pc_health_model$data[["y"]]),
                                         yrep = yrep_pc_health_model,
                                         group = pc_health_model$data[["wave"]],
                                         prob = 0.95,
                                         size = 0.5,
                                         fatten = 1,
                                         freq = FALSE) + # easier to compare groups (i.e., wave) with proportions, instead of counts
  ylim(c(0,1)) + theme(legend.position = "null") + ggtitle("Prior predictive check") +
  scale_x_continuous(name="Response options")
set.seed(NULL)

# posterior predictive check
set.seed(2021)

yrep_health_model <- posterior_predict(health_model)

postplot_health_model <- ppc_bars_grouped(y = as.numeric(health_model$data[["y"]]),
                                        yrep = yrep_health_model,
                                        group = health_model$data[["wave"]],
                                        prob = 0.95,
                                        size = 0.5,
                                        fatten = 1,
                                        freq = FALSE) + # easier to compare groups (i.e., wave) with proportions, instead of counts
  ylim(c(0,1)) + theme(legend.position = "null") + ggtitle("Posterior predictive check") +
  scale_x_continuous(name="Response options")
set.seed(NULL)

# plot
cairo_pdf("health_model_ppc.pdf", width = 7, height = 5.5) # start print to pdf

(priorplot_health_model + postplot_health_model + patchwork::plot_layout(ncol = 1, nrow = 2))

dev.off() # end print to pdf

## m_mrp:

# prior predictive check
set.seed(2021)

yrep_pc_m_mrp <- posterior_predict(priorcheck_m_mrp, resp = "y")

priorplot_m_mrp <- ppc_bars_grouped(y = as.numeric(priorcheck_m_mrp$data[["y"]]),
                                     yrep = yrep_pc_m_mrp,
                                     group = priorcheck_m_mrp$data[["wave"]],
                                     prob = 0.95,
                                     size = 0.5,
                                     fatten = 1,
                                     freq = FALSE) + # easier to compare groups (i.e., wave) with proportions, instead of counts
  ylim(c(0,1)) + theme(legend.position = "null") + ggtitle("Prior predictive check") +
  scale_x_continuous(name="Response options")
set.seed(NULL)

# posterior predictive check
set.seed(2021)

yrep_m_mrp <- posterior_predict(m_mrp, resp = "y")

postplot_m_mrp <- ppc_bars_grouped(y = as.numeric(m_mrp$data[["y"]]),
                                    yrep = yrep_m_mrp,
                                    group = m_mrp$data[["wave"]],
                                    prob = 0.95,
                                    size = 0.5,
                                    fatten = 1,
                                    freq = FALSE) + # easier to compare groups (i.e., wave) with proportions, instead of counts
  ylim(c(0,1)) + theme(legend.position = "null") + ggtitle("Posterior predictive check") +
  scale_x_continuous(name="Response options")
set.seed(NULL)

# plot
cairo_pdf("m_mrp_ppc.pdf", width = 7, height = 5.5) # start print to pdf

(priorplot_m_mrp + postplot_m_mrp + patchwork::plot_layout(ncol = 1, nrow = 2))

dev.off() # end print to pdf

## m1_imp:

# prior predictive check
set.seed(2021)

yrep_pc_m1_imp <- posterior_predict(m1_imp_ppc, resp = "y")

priorplot_m1_imp <- ppc_bars_grouped(y = as.numeric(m1_imp_ppc$data[["y"]]),
                                 yrep = yrep_pc_m1_imp,
                                 group = m1_imp_ppc$data[["wave"]],
                                 prob = 0.95,
                                 size = 0.5,
                                 fatten = 1,
                                 freq = FALSE) + # easier to compare groups (i.e., wave) with proportions, instead of counts
  ylim(c(0,1)) + theme(legend.position = "null") + ggtitle("Prior predictive check") +
  scale_x_continuous(name="Response options")
set.seed(NULL)

# posterior predictive check
set.seed(2021)

yrep_m1_imp <- posterior_predict(m1_imp, resp = "y")

postplot_m1_imp <- ppc_bars_grouped(y = as.numeric(m1_imp$data[["y"]]),
                                yrep = yrep_m1_imp,
                                group = m1_imp$data[["wave"]],
                                prob = 0.95,
                                size = 0.5,
                                fatten = 1,
                                freq = FALSE) + # easier to compare groups (i.e., wave) with proportions, instead of counts
  ylim(c(0,1)) + theme(legend.position = "null") + ggtitle("Posterior predictive check") +
  scale_x_continuous(name="Response options")
set.seed(NULL)

# plot
cairo_pdf("m1_imp_ppc.pdf", width = 7, height = 5.5) # start print to pdf

(priorplot_m1_imp + postplot_m1_imp + patchwork::plot_layout(ncol = 1, nrow = 2))

dev.off() # end print to pdf

### POSTERIOR PREDICTED CUMULATIVE PROBABILITIES

## m0:
set.seed(2021)

m0_fitted <- m0$data %>%
  modelr::data_grid(wave = modelr::seq_range(1:3, by = 1)) %>%
  add_epred_draws(m0, ndraws = 300, newdata = ., re_formula = NA, robust = TRUE) %>%
  mutate(category = factor(as.numeric(.category))) 

set.seed(NULL)

m0_fitted <- m0_fitted %>%
  group_by(.draw, category) %>%
  mutate(indices = cur_group_id()) %>%
  ungroup()

m0_fitted$c <- 0 # for facet_wrap() later

## m0_complete:

set.seed(2021)

m0_complete_fitted <- m0_complete$data %>%
  modelr::data_grid(wave = modelr::seq_range(1:3, by = 1)) %>%
  add_epred_draws(m0_complete, ndraws = 300, newdata = ., re_formula = NA, robust = TRUE) %>%
  mutate(category = factor(as.numeric(.category))) 

set.seed(NULL)

m0_complete_fitted <- m0_complete_fitted %>%
  group_by(.draw, category) %>%
  mutate(indices = cur_group_id()) %>%
  ungroup()

m0_complete_fitted$c <- 1 # for facet_wrap() later

m0_plot_df <- rbind(m0_fitted, m0_complete_fitted)

cairo_pdf("m0_cumprob.pdf", width = 8.5, height = 4.5) # start print to pdf

post_est_plot <- m0_plot_df %>%
  ggplot(aes(x = wave, 
             y = .epred, 
             color = .category, 
             group = indices))  + 
  geom_line(alpha = .1) +
  viridis::scale_color_viridis(discrete = TRUE, 
                               option = "D",
                               alpha = 1,
                               name="How important is \nreligion in your life?",
                               labels=c("Not at all important", "Not very important", "Somewhat important", "Very important")) +
  facet_wrap(~ factor(c, levels = c("0", "1"))) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(size=12, hjust = 0.5),
        axis.text = element_text( size = 10),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        plot.margin = margin(1, 1, 1, 1)) + 
  scale_x_continuous(name="Wave", breaks=c(1,2,3), labels= c(1,2,3)) +
  scale_y_continuous(name="Probability", breaks=c(0, .2, .4, .6, .8, 1), limits = c(0,1)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + # controls the legend alpha separately from the geom_line() alpha
  ggtitle("Posterior predictions from unadjusted models: listwise deletion (left) vs. complete-cases (right)")

post_est_plot

dev.off() # end to print pdf

## m1_complete

set.seed(2021)

m1_complete_fitted <- m1_complete$data %>%
  modelr::data_grid(wave = modelr::seq_range(1:3, by = 1), age.c = 0, gender = c(2,1), edu = 4, health = 5, hinc = 6) %>%
  add_epred_draws(m1_complete, ndraws = 300, newdata = ., re_formula = NA, robust = TRUE) %>%
  mutate(category = factor(as.numeric(.category))) 

set.seed(NULL)

m1_complete_fitted <- m1_complete_fitted %>%
  group_by(.draw, category) %>%
  mutate(indices = cur_group_id()) %>%
  ungroup()

cairo_pdf("m1_complete_cumprob.pdf", width = 8.5, height = 4.5) # start print to pdf

m1_complete_fitted %>%
  ggplot(aes(x = wave, 
             y = .epred, 
             color = category, 
             group = indices))  + 
  facet_wrap(~factor(gender, levels = c("2","1")), nrow = 1) + 
  geom_line(alpha = .1) +
  viridis::scale_color_viridis(discrete = TRUE, 
                               option = "D",
                               alpha = 1,
                               name="How important is \nreligion in your life?",
                               labels=c("Not at all important", "Not very important", "Somewhat important", "Very important")) + 
  theme_bw() +
  theme(
    panel.border = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    plot.title = element_text(size=12, hjust = 0.5),
    axis.text = element_text( size = 10),
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    plot.margin = margin(1, 1, 1, 1)) + # top, right, bottom, left
  scale_x_continuous(name="Wave", breaks=c(1,2,3)) +
  scale_y_continuous(name="Probability", breaks=c(0, .2, .4, .6, .8, 1), limits = c(0,1)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  ggtitle("Posterior predictions from adjusted complete-cases model: men (left) and women (right)")

dev.off() # end to print pdf

## m1_imp

set.seed(2021)

m1_imp_fitted <- m1_imp$data %>%
  modelr::data_grid(wave = modelr::seq_range(1:3, by = 1), age.c = 0, gender = c(2,1), edu = 4, health = 5, hinc = 6) %>%
  add_epred_draws(m1_imp, ndraws = 300, newdata = ., re_formula = NA, robust = TRUE, resp = "y") %>%
  mutate(category = factor(as.numeric(.category))) 

set.seed(NULL)

m1_imp_fitted <- m1_imp_fitted %>%
  group_by(.draw, category) %>%
  mutate(indices = cur_group_id()) %>%
  ungroup()

cairo_pdf("m1_imp_cumprob.pdf", width = 8.5, height = 4.5) # start print to pdf

m1_imp_fitted %>%
  ggplot(aes(x = wave, 
             y = .epred, 
             color = category, 
             group = indices))  + 
  facet_wrap(~factor(gender, levels = c("2","1")), nrow = 1) + 
  geom_line(alpha = .1) +
  viridis::scale_color_viridis(discrete = TRUE, 
                               option = "D",
                               alpha = 1,
                               name="How important is \nreligion in your life?",
                               labels=c("Not at all important", "Not very important", "Somewhat important", "Very important")) + 
  theme_bw() +
  theme(
    panel.border = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    plot.title = element_text(size=12, hjust = 0.5),
    axis.text = element_text( size = 10),
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    plot.margin = margin(1, 1, 1, 1)) + # top, right, bottom, left
  scale_x_continuous(name="Wave", breaks=c(1,2,3)) +
  scale_y_continuous(name="Probability", breaks=c(0, .2, .4, .6, .8, 1), limits = c(0,1)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  ggtitle("Posterior predictions from imputation model: men (left) and women (right)")

dev.off() #end print to pdf

## m_mrp

# credit to:
# https://www.monicaalexander.com/posts/2019-08-07-mrp/
# https://rohanalexander.com/posts/2019-12-04-getting_started_with_mrp/
# https://osf.io/preprints/socarxiv/3v5g7/

# construct external dataset to use for poststratification:
# https://appsso.eurostat.ec.europa.eu/nui/submitViewTableAction.do
# https://appsso.eurostat.ec.europa.eu/nui/show.do?dataset=lfsa_pgaed

# load external data for poststratification
eurostat <- read.csv("eurostat.csv", sep = ";")

# create data for poststratification
post_df <- with(eurostat, data.frame(wave = rep(seq(1, 3, by = 1), each = 18), 
                                     gender = rep(SEX, 3), 
                                     age = rep(AGE, 3), 
                                     edu = rep(ISCED11, 3), 
                                     value = rep(Value, 3))) # number of people with given combination of demographic characteristics

# bin demographic variables
post_df$gender[post_df$gender == "Females"] <- 1 # women
post_df$gender[post_df$gender == "Males"] <- 2 # men

post_df$age[post_df$age == "Y15-24"] <- 1 # 15-24
post_df$age[post_df$age == "Y25-49"] <- 2 # 25-49
post_df$age[post_df$age == "Y50-74"] <- 3 # 50-74

post_df$edu[post_df$edu == "Less than primary, primary and lower secondary education (levels 0-2)"] <- 1 # primary to lower secondary
post_df$edu[post_df$edu == "Upper secondary and post-secondary non-tertiary education (levels 3 and 4)"] <- 2 # # upper secondary and post-secondary non-tertiary education
post_df$edu[post_df$edu == "Tertiary education (levels 5-8)"] <- 3 # tertiary education

str(post_df)

# Table S1: table of external poststrat data

eurostat$ISCED11[eurostat$ISCED11 == "Less than primary, primary and lower secondary education (levels 0-2)"] <- "Levels 0-2" # primary to lower secondary
eurostat$ISCED11[eurostat$ISCED11 == "Upper secondary and post-secondary non-tertiary education (levels 3 and 4)"] <- "Levels 3-4" # # upper secondary and post-secondary non-tertiary education
eurostat$ISCED11[eurostat$ISCED11 == "Tertiary education (levels 5-8)"] <- "Levels 5-8" # tertiary education

eurostat$SEX[eurostat$SEX == "Females"] <- "Female" # women
eurostat$SEX[eurostat$SEX == "Males"] <- "Male" # men

post_tab <- eurostat  %>%
  group_by(SEX, AGE, ISCED11)  %>% 
  summarize(N = round(Value*1000)) %>% 
  distinct() %>%
  ungroup() 

colnames(post_tab) <- c("Gender", "Age", "Education", "N")

print(xtable(post_tab, 
             align = c("l", "l", "c", "c", "c"), digits = c(rep(0, 5)),
             caption = "Poststratification matrix"), include.rownames = FALSE)

# overall poststratification
ps_prop_pop <- post_df  %>% # get proportions for each combination of demographic variables
  group_by(wave) %>%
  mutate(prop = value/sum(value)) %>%
  ungroup()

ps_prop_pop

set.seed(2021)

post_est_pop_fitted <- m_mrp %>%
  add_epred_draws(newdata=ps_prop_pop, ndraws = 300, re_formula = ~ (1 | age) + (1 | gender) + (1 | edu), allow_new_levels = TRUE) %>%
  rename(estimate = .epred) %>% 
  mutate(estimate_prop = estimate*prop) %>% 
  group_by(wave, .draw, .category) %>% 
  summarise(estimate_sum = sum(estimate_prop)) %>% 
  group_by(.draw, .category) %>%
  mutate(indices = cur_group_id()) %>%
  ungroup() %>%
  rename(.epred = estimate_sum)

set.seed(NULL)

post_est_pop_fitted

post_est_pop_fitted$ps <- 1 # for plotting

# compare with unadjusted model (m0)
set.seed(2021)

no_post_est_pop_fitted <- m0 %>%
  add_epred_draws(newdata = data.frame(wave=seq(1, 3, by = 1)), re_formula = NA, ndraws = 300, allow_new_levels = TRUE) %>% 
  group_by(wave, .draw, .category) %>% 
  summarise(.epred = .epred) %>% 
  group_by(.draw, .category) %>%
  mutate(indices = cur_group_id()) %>%
  ungroup()

set.seed(NULL)

no_post_est_pop_fitted

no_post_est_pop_fitted$ps <- 0 # for plotting

# plot poststratified model predictions vs. unadjusted model (m0)
post_est <- rbind(post_est_pop_fitted, no_post_est_pop_fitted)

cairo_pdf("mrp_cumprob.pdf", width = 8.5, height = 4.5) # start print to pdf

post_est_plot <- post_est %>%
  ggplot(aes(x = wave, 
             y = .epred, 
             color = .category, 
             group = indices))  + 
  geom_line(alpha = .1) +
  viridis::scale_color_viridis(discrete = TRUE, 
                               option = "D",
                               alpha = 1,
                               name="How important is \nreligion in your life?",
                               labels=c("Not at all important", "Not very important", "Somewhat important", "Very important")) +
  facet_wrap(~ factor(ps, levels = c("1", "0"))) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(size=12, hjust = 0.5),
        axis.text = element_text( size = 10),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        plot.margin = margin(1, 1, 1, 1)) + 
  scale_x_continuous(name="Wave", breaks=c(1,2,3), labels= c(1,2,3)) +
  scale_y_continuous(name="Probability", breaks=c(0, .2, .4, .6, .8, 1), limits = c(0,1)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + # controls the legend alpha separately from the geom_line() alpha
  ggtitle("Posterior predictions with (left) and without (right) poststratification")

post_est_plot

dev.off() # end print to pdf

### INTERACTION PLOTS

## hinc_model:

hinc_c <-  as.data.frame(modelr::data_grid(d_full, wave = c(1,2,3), hinc = c(1,2,3,4,5,6,7,8,9,10,11))) # combination of covariates to condition on
hinc_cond <- conditional_effects(hinc_model, effect = "wave", conditions = hinc_c, # effect goes on the x-axis
                            robust = TRUE, # TRUE = posterior medians (instead of means)
                            categorical = TRUE, # set categorical = F to see (improper) continuous interaction
                            prob = .95, # toggle interval width
                            re_formula = NA) # include (NULL) or not (NA) random effects

hinc_labels <- c("1" = "Less than 100.000 kr.", "2" = "100.000 to 199.999 kr.", 
                "3" = "200.000 to 299.999 kr.", "4" = "300.000 to 399.999 kr.",
                "5" = "400.000 to 499.999 kr.", "6" = "500.000 to 599.999 kr.",
                "7" = "600.000 to 699.999 kr.", "8" = "700.000 to 799.999 kr.",
                "9" = "800.000 to 899.999 kr.", "10" = "900.000 to 999.999 kr.",
                "11" = "1.000.000 kr. or more")

cairo_pdf("hinc_cumprob.pdf", width = 8.5, height = 8.5) # start print to pdf

plot(hinc_cond, plot = FALSE)[[1]] +
  viridis::scale_color_viridis(discrete = TRUE, 
                               option = "D",
                               alpha = .8,
                               name="How important is \nreligion in your life?",
                               labels=c("Not at all important", "Not very important", "Somewhat important", "Very important")) + 
  viridis::scale_fill_viridis(discrete = TRUE, 
                              option = "D",
                              alpha = .8,
                              name="How important is \nreligion in your life?",
                              labels=c("Not at all important", "Not very important", "Somewhat important", "Very important")) +
  facet_wrap(~ factor(hinc, levels = c(1,2,3,4,5,6,7,8,9,10,11)), labeller = as_labeller(hinc_labels)) +
  theme_minimal() +
  theme(legend.position = c(0.9, 0.15)) +
  scale_x_continuous(breaks = c(1,2,3)) +
  scale_y_continuous(name="Probability", breaks=c(0, .2, .4, .6, .8, 1), limits = c(0,1)) + 
  ggtitle("Annual household income")

dev.off()

## health_model:

health_c <-  as.data.frame(modelr::data_grid(d_full, wave = c(1,2,3), health = c(1,2,3,4,5))) # combination of covariates to condition on
health_cond <- conditional_effects(health_model, effect = "wave", conditions = health_c, # effect goes on the x-axis
                                 robust = TRUE, # TRUE = posterior medians (instead of means)
                                 categorical = TRUE, # set categorical = F to see (improper) continuous interaction
                                 prob = .95, # toggle interval width
                                 re_formula = NA) # include (NULL) or not (NA) random effects

health_labels <- c("1" = "''Very bad''", "2" = "''Bad''", 
                 "3" = "''Okay''", "4" = "''Good''",
                 "5" = "''Very good''")

cairo_pdf("health_cumprob.pdf", width = 8.5, height = 8.5) # start print to pdf

plot(health_cond, plot = FALSE)[[1]] +
  viridis::scale_color_viridis(discrete = TRUE, 
                               option = "D",
                               alpha = .8,
                               name="How important is \nreligion in your life?",
                               labels=c("Not at all important", "Not very important", "Somewhat important", "Very important")) + 
  viridis::scale_fill_viridis(discrete = TRUE, 
                              option = "D",
                              alpha = .8,
                              name="How important is \nreligion in your life?",
                              labels=c("Not at all important", "Not very important", "Somewhat important", "Very important")) +
  facet_wrap(~ factor(health, levels = c(1,2,3,4,5)), labeller = as_labeller(health_labels)) +
  theme_minimal() +
  theme(legend.position = c(0.85, 0.25)) +
  scale_x_continuous(breaks = c(1,2,3)) +
  scale_y_continuous(name="Probability", breaks=c(0, .2, .4, .6, .8, 1), limits = c(0,1)) + 
  ggtitle("''How would you describe your current health?''")

dev.off()

## exposure models

exp_labels <- c("1" = "''No, but have not been tested''", "2" = "''No, and have tested negative''", 
                "3" = "''Yes, have tested positive but was not very ill''", "4" = "''Yes, have tested positive and was very ill''",
                "5" = "''Yes, was hospitalized''")

# m_self:
self_c <-  as.data.frame(modelr::data_grid(df_exp, wave = c(1,2,3), exp_self = c(1,2,3,4,5))) # combination of covariates to condition on
self_cond <- conditional_effects(m_self, effect = "wave", conditions = self_c, # effect goes on the x-axis
                                      robust = TRUE, # TRUE = posterior medians (instead of means)
                                      categorical = TRUE, # set categorical = F to see (improper) continuous interaction
                                      prob = .95, # toggle interval width
                                      re_formula = NA) # include (NULL) or not (NA) random effects

cairo_pdf("self_cumprob.pdf", width = 8.5, height = 8.5) # start print to pdf

plot(self_cond, plot = FALSE)[[1]] +
  viridis::scale_color_viridis(discrete = TRUE, 
                               option = "D",
                               alpha = .8,
                               name="How important is \nreligion in your life?",
                               labels=c("Not at all important", "Not very important", "Somewhat important", "Very important")) + 
  viridis::scale_fill_viridis(discrete = TRUE, 
                              option = "D",
                              alpha = .8,
                              name="How important is \nreligion in your life?",
                              labels=c("Not at all important", "Not very important", "Somewhat important", "Very important")) +
  facet_wrap("exp_self", labeller = as_labeller(exp_labels)) +
  theme_minimal() +
  theme(legend.position = c(0.85, 0.225)) +
  scale_x_continuous(breaks = c(1,2,3)) +
  scale_y_continuous(name="Probability", breaks=c(0, .2, .4, .6, .8, 1), limits = c(0,1)) + 
  ggtitle("''I've been ill [with corona virus] myself''")

dev.off()

# m_household:
household_c <-  as.data.frame(modelr::data_grid(df_exp, wave = c(1,2,3), exp_household = c(1,2,3,4,5))) # combination of covariates to condition on
household_cond <- conditional_effects(m_household, effect = "wave", conditions = household_c, # effect goes on the x-axis
                                 robust = TRUE, # TRUE = posterior medians (instead of means)
                                 categorical = TRUE, # set categorical = F to see (improper) continuous interaction
                                 prob = .95, # toggle interval width
                                 re_formula = NA) # include (NULL) or not (NA) random effects

cairo_pdf("household_cumprob.pdf", width = 8.5, height = 8.5) # start print to pdf

plot(household_cond, plot = FALSE)[[1]] +
  viridis::scale_color_viridis(discrete = TRUE, 
                               option = "D",
                               alpha = .8,
                               name="How important is \nreligion in your life?",
                               labels=c("Not at all important", "Not very important", "Somewhat important", "Very important")) + 
  viridis::scale_fill_viridis(discrete = TRUE, 
                              option = "D",
                              alpha = .8,
                              name="How important is \nreligion in your life?",
                              labels=c("Not at all important", "Not very important", "Somewhat important", "Very important")) +
  facet_wrap("exp_household", labeller = as_labeller(exp_labels)) +
  theme_minimal() +
  theme(legend.position = c(0.85, 0.225)) +
  scale_x_continuous(breaks = c(1,2,3)) +
  scale_y_continuous(name="Probability", breaks=c(0, .2, .4, .6, .8, 1), limits = c(0,1)) + 
  ggtitle("''Someone in my household is/has been ill [with corona virus]''")

dev.off()

# m_family:
family_c <-  as.data.frame(modelr::data_grid(df_exp, wave = c(1,2,3), exp_family = c(1,2,3,4,5))) # combination of covariates to condition on
family_cond <- conditional_effects(m_family, effect = "wave", conditions = family_c, # effect goes on the x-axis
                                      robust = TRUE, # TRUE = posterior medians (instead of means)
                                      categorical = TRUE, # set categorical = F to see (improper) continuous interaction
                                      prob = .95, # toggle interval width
                                      re_formula = NA) # include (NULL) or not (NA) random effects

cairo_pdf("family_cumprob.pdf", width = 8.5, height = 8.5) # start print to pdf

plot(family_cond, plot = FALSE)[[1]] +
  viridis::scale_color_viridis(discrete = TRUE, 
                               option = "D",
                               alpha = .8,
                               name="How important is \nreligion in your life?",
                               labels=c("Not at all important", "Not very important", "Somewhat important", "Very important")) + 
  viridis::scale_fill_viridis(discrete = TRUE, 
                              option = "D",
                              alpha = .8,
                              name="How important is \nreligion in your life?",
                              labels=c("Not at all important", "Not very important", "Somewhat important", "Very important")) +
  facet_wrap("exp_family", labeller = as_labeller(exp_labels)) +
  theme_minimal() +
  theme(legend.position = c(0.85, 0.225)) +
  scale_x_continuous(breaks = c(1,2,3)) +
  scale_y_continuous(name="Probability", breaks=c(0, .2, .4, .6, .8, 1), limits = c(0,1)) + 
  ggtitle("''Someone in my family is/has been ill [with corona virus]''")

dev.off()

# m_relations:

relations_c <-  as.data.frame(modelr::data_grid(df_exp, wave = c(1,2,3), exp_relations = c(1,2,3,4,5))) # combination of covariates to condition on
relations_cond <- conditional_effects(m_relations, effect = "wave", conditions = relations_c, # effect goes on the x-axis
                            robust = TRUE, # TRUE = posterior medians (instead of means)
                            categorical = TRUE, # set categorical = F to see (improper) continuous interaction
                            prob = .95, # toggle interval width
                            re_formula = NA) # include (NULL) or not (NA) random effects

cairo_pdf("relations_cumprob.pdf", width = 8.5, height = 8.5) # start print to pdf

plot(relations_cond, plot = FALSE)[[1]] +
  viridis::scale_color_viridis(discrete = TRUE, 
                               option = "D",
                               alpha = .8,
                               name="How important is \nreligion in your life?",
                               labels=c("Not at all important", "Not very important", "Somewhat important", "Very important")) + 
  viridis::scale_fill_viridis(discrete = TRUE, 
                              option = "D",
                              alpha = .8,
                              name="How important is \nreligion in your life?",
                              labels=c("Not at all important", "Not very important", "Somewhat important", "Very important")) +
  facet_wrap("exp_relations", labeller = as_labeller(exp_labels)) +
  theme_minimal() +
  theme(legend.position = c(0.85, 0.225)) +
  scale_x_continuous(breaks = c(1,2,3)) +
  scale_y_continuous(name="Probability", breaks=c(0, .2, .4, .6, .8, 1), limits = c(0,1)) + 
  ggtitle("''Someone in my relations is/has been ill [with corona virus]''")

dev.off()
