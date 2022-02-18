###########################
### MAIN TABLES & PLOTS ###
###########################

### Install packages, their dependencies, and specific versions locally to ensure a reproducible package environment
renv::use(
  "BH@1.75.0-0",
  "Brobdingnag@1.2-6",
  "DBI@1.1.2",
  "DT@0.20",
  "HDInterval@0.2.2",
  "JointAI@1.0.3",
  "MASS@7.3-54",
  "Matrix@1.3-4",
  "R6@2.5.1",
  "RColorBrewer@1.1-2",
  "Rcpp@1.0.8",
  "RcppArmadillo@0.10.7.5.0",
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
  "brio@1.1.2",
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
  "diffobj@0.3.5",
  "digest@0.6.29",
  "distributional@0.3.0",
  "dplyr@1.0.7",
  "dtplyr@1.2.1",
  "dygraphs@1.1.1.6",
  "ellipse@0.4.2",
  "ellipsis@0.3.2",
  "evaluate@0.14",
  "fansi@1.0.2",
  "farver@2.1.0",
  "fastmap@1.1.0",
  "fftwtools@0.9-11",
  "fontawesome@0.2.2",
  "forcats@0.5.1",
  "foreach@1.5.1",
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
  "iterators@1.0.13",
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
  "mathjaxr@1.4-0",
  "matrixStats@0.61.0",
  "mcmcse@1.5-0",
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
  "pillar@1.7.0",
  "pkgbuild@1.3.1",
  "pkgconfig@2.0.3",
  "pkgload@1.2.4",
  "plyr@1.8.6",
  "posterior@1.2.0",
  "praise@1.0.0",
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
  "rjags@4-12",
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
  "survival@3.2-13",
  "svUnit@1.0.6",
  "sys@3.4",
  "tensorA@0.36.2",
  "testthat@3.1.2",
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
  "waldo@0.3.1",
  "withr@2.4.3",
  "xfun@0.29",
  "xml2@1.3.2",
  "xtable@1.8-4",
  "xts@0.12.1",
  "yaml@2.2.2",
  "zoo@1.8-9"
)

### Plotting 
library(modelr) # for plotting
library(xtable) # for printing TeX tables
library(brms) # for implementing Bayesian multilevel analysis
library(tidyverse) # wrangling, plotting, etc.
library(tidybayes) # for preparing and visualizing prior and posterior draws
library(bayesplot) # for prior and posterior predictive plots
library(JointAI) # for overview of missing values

# renv::embed() # to initialize and embed lock file in script; generates the renv::use(...) above

# load all fitted models
if (file.exists("pandemic_panel_fits.RData")) base::load(file = "pandemic_panel_fits.RData")

##################################
### Overview of missing values ###
##################################

JointAI::plot_all(d_full) # overall
JointAI::plot_all(d_full[d_full$wave==1,]) # wave 1
JointAI::plot_all(d_full[d_full$wave==2,]) # wave 2
JointAI::plot_all(d_full[d_full$wave==3,]) # wave 3

####################
### Demographics ###
####################

### TABLE 1: demographic table

(demo.perc <- d_full[complete.cases(d_full),] %>%
   group_by(wave) %>%
   summarize(
     N = length(unique(id)),
     mAge = mean(age),
     sdAge = sd(age),
     "Gender (women)" = sum(as.numeric(gender[gender==1]))/length(unique(id))*100,
     Education = t(table(as.numeric(edu)))/length(unique(id))*100,
     "Household income" = t(table(hinc))/length(unique(id))*100
   ) %>% 
   t() %>%
   as.data.frame()
)

demo.perc <- demo.perc[-1,]

digmat.perc <- matrix(c(0,0,0,0, rep(c(0, 1, 1, 1),(20))), 
                      nrow = 21, ncol=4, byrow=TRUE)

xtable(demo.perc, digits = digmat.perc, align = c("l", "c", "c", "c"),
       caption = "Demographic overview for each wave. Mean age and standard deviation in parantheses. Percentages for gender, educational and household income levels. Complete-case sample.")

length(unique(d_full[complete.cases(d_full),]$id)) # total number of unique id's, identical to data used for fitting m1

#################################
### PLOTS FOR MAIN MANUSCRIPT ###
#################################

### FIGURE 2: Ordered id coef plot
df_ordplot <- as.data.frame(coef(m1, summary = TRUE, pars = "mowave", 
                                 probs = c(0.25, 0.75, 0.025, 0.975))[["id"]]) %>% # include both 50% and 95% intervals
  rownames_to_column("id") %>% 
  as_tibble() %>% 
  rename(wave_slope = Estimate.mowave) %>% 
  as.data.frame()

df_ordplot <- df_ordplot[order(df_ordplot$wave_slope),] # order id's for their beta coef.
df_ordplot$order_id <- 1:nrow(df_ordplot) # give each participant a new id, ordered for their beta coef.

cairo_pdf("Figure2.pdf",
          width = 14, height = 8) # start print to pdf

ggplot(df_ordplot, aes(x=reorder(order_id, wave_slope), y=wave_slope)) + 
  geom_linerange(aes(ymin = Q2.5.mowave, ymax = Q97.5.mowave, color = wave_slope), alpha = .60) + # 95% CI intervals
  geom_point(size = .1) + # posterior means as points
  geom_linerange(aes(ymin = Q25.mowave, ymax = Q75.mowave), alpha = .2) + # 50% CI intervals
  geom_hline(yintercept = 0, linetype = 'dotted', size = .5) +
  viridis::scale_color_viridis(discrete = FALSE, 
                               option = "D",
                               alpha = 1) + 
  viridis::scale_fill_viridis(discrete = FALSE, 
                              option = "D",
                              alpha = 1) +
  scale_y_continuous(name="Slope coefficient (wave)", limits = c(-4,4)) +
  scale_x_discrete(name="Participant number", 
                   breaks=c(1,max(df_ordplot$order_id))) +
  theme_minimal() +
  theme(axis.text = element_text( size = 25),
        legend.position = "none",
        axis.title.x = element_text(size=30),
        axis.title.y = element_text(size=30),
        plot.margin = unit(c(.1,1,.1,.25), "cm"))

dev.off() # end print to pdf

### FIGURE 3: Posterior predictive draws

# credit to: https://octavio.me/posts/2021-07-23-ordinal-viz/

set.seed(2021)

m1_fitted <- d_full %>%
  modelr::data_grid(wave = seq_range(1:3, by = 1), age.c = 0, gender = c(2,1), edu = 4, health = 5, hinc = 6) %>%
  add_epred_draws(m1, ndraws = 300, newdata = ., re_formula = NA, robust = TRUE, resp = "y") %>%
  mutate(category = factor(as.numeric(.category))) 

set.seed(NULL)

m1_fitted <- m1_fitted %>%
  group_by(.draw, category) %>%
  mutate(indices = cur_group_id()) %>%
  ungroup()

cairo_pdf("Figure3.pdf",
          width = 8.5, height = 4.5) # start print to pdf

m1_fitted %>%
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
  ggtitle("Posterior predictions for men (left) and women (right)") 

dev.off() # end print to pdf

### FIGURE 4: Posterior predictions for 15 individuals

id_plot_df <- d_full # copy full dataset

# only complete-case individuals so that the model can predict for covariates
id_plot_df_complete <- d_full[complete.cases(id_plot_df[4:9]),] %>% group_by(id) %>% filter(n()>2) %>% as.data.frame()

set.seed(2021) # make random sampling reproducible -- change seed to sample other IDs
id_plot_id <- sample(unique(id_plot_df_complete$id), 15, replace = FALSE) # sample randomly 15 id numbers
set.seed(NULL)

id_plot_df_select <- id_plot_df_complete %>%
  filter(id %in% id_plot_id) %>% # select the 15 participants
  mutate(y = as.numeric(y)) %>% 
  replace_na(list(y = 9)) # replace NA's with impossible value, so that rows with NAs are not excluded, allowing us to illustrate predictions for unobserved values

cairo_pdf("Figure4.pdf", width = 6, height = 4.5) # start print to pdf

set.seed(2021)

id_plot_ribbon_m1 <- id_plot_df_select %>%
  dplyr::group_by(id) %>%
  # set predictions for each participant's combination of covariates
  with(data.frame(id = id, wave = wave,
                  health = health, gender = gender, edu = edu, 
                  hinc = hinc, age.c = age.c)) %>%
  add_predicted_draws(newdata = ., object=m1, ndraws = 200, re_formula = NULL) %>%
  ggplot() +
  stat_lineribbon(aes(x = wave, y = as.numeric(.prediction), group = paste(id)), .width = c(.5, .95), alpha = .75, color = "#08519C",
                  point_interval = median_qi) + # could also plot the highest posterior density, using: mode_hdci
  geom_point(data = id_plot_df_select, aes(x = wave, y = y, group = paste(id)), color = "black",
             alpha = .8, size = 2) +
  facet_wrap(~id, nrow = 3) +
  scale_fill_brewer() +
  scale_x_continuous(name="Wave", breaks=c(1,2,3), labels = c(1,2,3)) +
  scale_y_continuous(name="Response options", breaks=c(1,2,3,4), limits = c(.8,4.2)) + 
  theme_minimal() +
  theme(legend.position = "blank",
        strip.background = element_blank(), # remove grey background for title
        axis.title = element_text( size = 11), # change size and bold font of global axis titles
        axis.text = element_text( size = 10), # axes text
        strip.text = element_text(size = 8), # ID text size
        panel.grid.minor = element_blank(), # remove minor grid lines
        panel.border = element_rect(colour = "darkgrey", fill=NA, size=.5)) # border line size and color

id_plot_ribbon_m1 # throws an error ("Removed 5 rows containing missing values (geom_point).") which can be ignored (it's just removing NAs)

set.seed(NULL)

dev.off() # end print to pdf

##################################
### Table 3: Coefficient table ###
##################################

coefnoadj <- c("Intercept[1]", "Intercept[2]", "Intercept[3]", "mowave", "mowave:mohinc", "mowave:mohealth")
coefadj <- c("Intercept[1]", "Intercept[2]", "Intercept[3]", "mowave", "mowave:mohinc", "mowave:mohealth")
coefimp <- c("y_Intercept[1]", "y_Intercept[2]", "y_Intercept[3]", "y_mowave", "y_mowave:mihinc", "y_mowave:mihealth")

fixef_m0 <- as.data.frame(round(fixef(m0, pars = coefnoadj), 2))
fixef_m1 <- as.data.frame(round(fixef(m1, pars = coefadj), 2))
fixef_m0_complete <- as.data.frame(round(fixef(m0_complete, pars = coefnoadj), 2))
fixef_m1_complete <- as.data.frame(round(fixef(m1_complete, pars = coefadj), 2))
fixef_health <- as.data.frame(round(fixef(health_model, pars = coefadj), 2))
fixef_hinc <- as.data.frame(round(fixef(hinc_model, pars = coefadj), 2))
fixef_imp <- as.data.frame(round(fixef(m1_imp, pars = coefimp), 2))

coef_list <- list(m0 = paste0(fixef_m0[,1], " [",fixef_m0[,3], ", ", fixef_m0[,4],"]"),
                            
                       m1 = paste0(fixef_m1[,1], " [",fixef_m1[,3], ", ", fixef_m1[,4],"]"),
                 
                       m0_complete = paste0(fixef_m0_complete[,1], " [",fixef_m0_complete[,3], ", ", fixef_m0_complete[,4],"]"),
                 
                       m1_complete = paste0(fixef_m1_complete[,1], " [",fixef_m1_complete[,3], ", ", fixef_m1_complete[,4],"]"),
                  
                       hinc_model = paste0(fixef_hinc[,1], " [",fixef_hinc[,3], ", ", fixef_hinc[,4],"]"),
                 
                       health_model = paste0(fixef_health[,1], " [",fixef_health[,3], ", ", fixef_health[,4],"]"),
                 
                       m1_imp = paste0(fixef_imp[,1], " [",fixef_imp[,3], ", ", fixef_imp[,4],"]"))

coef_list[[1]] <- append(coef_list[[1]], "--", after = 4)
coef_list[[1]] <- append(coef_list[[1]], "--", after = 5)
coef_list[[3]] <- append(coef_list[[3]], "--", after = 4)
coef_list[[3]] <- append(coef_list[[3]], "--", after = 5)
coef_list[[5]] <- append(coef_list[[5]], "--", after = 5)
coef_list[[6]] <- append(coef_list[[6]], "--", after = 4)

coef_tab <- as.data.frame(coef_list)

rownames(coef_tab) <- c("First threshold", "Second threshold", "Third threshold", "Wave", "Wave x H. income", "Wave x Health")

print(xtable(t(coef_tab), 
             align = c("l", "l", "l", "l", "l", "l", "l"),
             caption = "Key coefficients across main models"))

### END 
