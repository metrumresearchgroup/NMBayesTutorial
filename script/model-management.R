# The model-management.R file is intended to be a scratchpad for doing things
# like defining, submitting, tagging, etc. your models. There is no need to keep
# a "record in code" of these activities because they can all be reconstructed
# later via functions like `run_log()`, as demonstrated in `model-summary.Rmd`
#
# The `Model Management Demo` (rendered from the `model-management-demo.Rmd`
# file) shows code for a range of these activities at different stages in the
# modeling process. It exists purely for reference; the intent is _not_ for you
# to replicate the full narrative.
# https://ghe.metrumrg.com/pages/example-projects/bbr-nonmem-poppk-foce/model-management-demo
#
# This script assumes you have already installed and set up bbi. For details
# on getting set up with bbr, see:
# https://metrumresearchgroup.github.io/bbr/articles/getting-started.html#setup


library(bbr)
library(tidyverse)

setwd(file.path(rprojroot::find_rstudio_root_file(),"script"))
source("functions-model.R")


# define model dir and load tags
MODEL_DIR <- "../model/pk"
TAGS <- yaml::read_yaml("tags.yaml")


############################################
# CHECK FOR BBI INSTALLATION AND CONFIG 
############################################

# The code below checks that you have everything configured to begin modeling
# with bbr. This code can be deleted once everything is working correctly.

# To check if you have bbi installed and configured, run `bbi_version()`
bbi_version()
# If this doesn't return a version number, see 
# https://metrumresearchgroup.github.io/bbr/articles/getting-started.html#installing-bbi
# for installation details.


# The first time you are modeling in a new directory, you will need to "intialize" bbi.
# The bbi_init() function will create a bbi.yaml, with the default settings, in the 
# specified directory. 
# 
# To check you have initialized bbi, try to read in the `bbi.yaml` file.
file.path(MODEL_DIR, "bbi.yaml") %>% yaml::read_yaml() %>% names()
# If this errors, run `bbi_init()`:

bbi_init(.dir = MODEL_DIR,            # the directory to create the bbi.yaml in
         .nonmem_dir = "/opt/NONMEM", # location of NONMEM installation
         .nonmem_version = "nm75")  # default NONMEM version to use

# Note this only needs to be done _once for each folder_ you are modeling in. Once the bbi.yaml exists, 
# you will not need to run `bbi_init()` again unless you want to create another one; for example if you 
# move to modeling in a different directory.
#
# For more details on the `bbi.yaml` file and its usage, see:
# https://metrumresearchgroup.github.io/bbr/articles/getting-started.html#bbi-yaml-configuration-file

####Demo Example#####
mod1000 <- new_model(file.path(MODEL_DIR, 1000)) %>%
  add_tags(c(
    TAGS$one_compartment_absorption,
    TAGS$eta_cl,
    TAGS$eta_ka,
    TAGS$eta_v,
    TAGS$cov_cl_wt_fixed,
    TAGS$cov_v2_wt_fixed,
    TAGS$cov_q_wt_fixed,
    TAGS$cov_v3_wt_fixed,
    TAGS$proportional_ruv,
    TAGS$bayes
  ))
#Model object to be used with bbr
mod1000 <- read_model(file.path(MODEL_DIR, 1000))

#Submitting the initial chain stub file to generate
#starting values.
submit_model(
  mod1000,
  .bbi_args = list(overwrite = TRUE),
  .mode = "local"
)

#Quick examination of the sampled values in 1000.chn
read.table(file.path(MODEL_DIR, "1000/1000.chn"), header = TRUE) %>% 
  select(ITERATION:`OMEGA.3.3.`) %>% 
  pivot_longer(cols = -ITERATION) %>% 
  pivot_wider(names_from = "ITERATION")

#This function take the chain stub file (1000.ctl) and writes out 
#four additional control streams (appended with a numeral to denote 
#the chain number) where each control stream pulls from 
#the chain file (1000.chn) to use the values as starting points for 
#the sampling. The four files are then submitted and their output is 
#kept in the same directory.
run_chains(MODEL_DIR, 1000, .bbi_args = list(
  overwrite = TRUE, parallel = TRUE, threads = 2),
  .mode = "local")


####Pediatrics Example####
#Prior predictive simulation control stream
mod2000pps <- new_model(file.path(MODEL_DIR, '2000pps')) %>%
  add_tags(c(
    TAGS$bayes
  ))

submit_model(
  mod2000pps,
  .bbi_args = list(overwrite = TRUE),
  .mode = "local"
)


mod2000 <- new_model(file.path(MODEL_DIR, 2000)) %>%
  add_tags(c(
    TAGS$bayes
  ))

mod2000 <- read_model(file.path(MODEL_DIR, 2000))

submit_model(
  mod2000,
  .bbi_args = list(overwrite = TRUE),
  .mode = "local"
)

read.table(file.path(MODEL_DIR, "2000/2000.chn"), header = TRUE) %>% 
  select(ITERATION:`OMEGA.2.2.`) %>% 
  pivot_longer(cols = -ITERATION) %>% 
  pivot_wider(names_from = "ITERATION")

run_chains(MODEL_DIR, 2000, .bbi_args = list(
  overwrite = TRUE, parallel = TRUE, threads = 2),
  .mode = "local")

#####Example Sensitivity Analysis for Pediatric Example######
#Five models will be run
#One with KA fixed (2001)
#One with KA prior inflated 10x (2002)
#One with KA prior variance 100x (2003)
#One with KA prior mean increased 50% (2004)
#One with KA prior mean decreased 50% (2005)

mod2001 <- copy_model_from(mod2000,2001) %>%
  update_run_number()
edit_model(mod2001)

mod2002 <- copy_model_from(mod2000,2002) %>%
  update_run_number()
edit_model(mod2002)

mod2003 <- copy_model_from(mod2000,2003) %>%
  update_run_number()
edit_model(mod2003)

mod2004 <- copy_model_from(mod2000,2004) %>%
  update_run_number()
edit_model(mod2004)

mod2005 <- copy_model_from(mod2000,2005) %>%
  update_run_number()
edit_model(mod2005)

submit_models(
  list(mod2001,mod2002, mod2003,mod2004,mod2005),
  .bbi_args = list(overwrite = TRUE),
  .mode = "local"
)

run_chains(MODEL_DIR, 2001, .bbi_args = list(
  overwrite = TRUE, parallel = TRUE, threads = 2),
  .mode = "local")

