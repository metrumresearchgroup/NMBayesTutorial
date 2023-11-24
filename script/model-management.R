library(bbr)
library(bbr.bayes)
library(tidyverse)
library(here)

# define model dir and load tags
MODEL_DIR <- here("model/pk")
TAGS <- yaml::read_yaml("script/tags.yaml")


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

# Start with FOCE model (suppose model has been developed using FOCE to this point)

mod1000_foce <- new_model(file.path(MODEL_DIR, "1000-FOCE")) %>%
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
    TAGS$foce
  ))
#Model object to be used with bbr
mod1000_foce <- read_model(file.path(MODEL_DIR, "1000-FOCE"))

# Run the FOCE model
submit_model(
  mod1000_foce,
  .bbi_args = list(overwrite = TRUE),
  .mode = "local"
)

# Copy FOCE model to Bayes template model
mod1000 <- copy_model_as_nmbayes(
  mod1000_foce,
  file.path(MODEL_DIR, 1000),
  .inherit_tags = TRUE
) %>% 
  replace_tag(TAGS$foce, TAGS$bayes)

# This sets up a new control stream (1000.ctl) with the model mostly unchanged
# apart from the following:

# ; TODO: This model was copied by bbr.bayes::copy_model_as_nmbayes().
# ;       nmbayes models require a METHOD=CHAIN estimation record and a
# ;       METHOD=BAYES or METHOD=NUTS estimation record. The records
# ;       below are meant as a starting point.  At the very least, you
# ;       need to adjust the number of iterations (see NITER option in
# ;       the second $EST block), but please review all options
# ;       carefully.
# ;
# ;       See ?bbr.bayes::bbr_nmbayes and the NONMEM docs for details.
# $EST METHOD=CHAIN FILE=1000.chn NSAMPLE=4 ISAMPLE=0 SEED=1
# CTYPE=0 IACCEPT=0.3 DF=10 DFS=0
# 
# $EST METHOD=NUTS SEED=1 NBURN=250 NITER=NNNN
# AUTO=2 CTYPE=0 OLKJDF=2 OVARF=1
# NUTS_DELTA=0.95 PRINT=10 MSFO=1000.msf RANMETHOD=P PARAFPRINT=10000
# BAYES_PHI_STORE=1

# As per these instructions, we edit 1000.ctl which will act as a template for 4
# chain runs. The main change will be the addition of priors.

#Model object to be used with bbr
mod1000 <- read_model(file.path(MODEL_DIR, 1000))

# Run the model
submit_model(
  mod1000,
  .bbi_args = list(overwrite = TRUE, threads = 2)
)

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
  .bbi_args = list(overwrite = TRUE, threads = 2)
)

#####Example Sensitivity Analysis for Pediatric Example######
#Five models will be run
#One with KA fixed (2001)
#One with KA prior inflated 10x (2002)
#One with KA prior variance 100x (2003)
#One with KA prior mean increased 50% (2004)
#One with KA prior mean decreased 50% (2005)

mod2001 <- copy_model_from(mod2000,2001) %>%
  update_model_id()
edit_model(mod2001)

mod2002 <- copy_model_from(mod2000,2002) %>%
  update_model_id()
edit_model(mod2002)

mod2003 <- copy_model_from(mod2000,2003) %>%
  update_model_id()
edit_model(mod2003)

mod2004 <- copy_model_from(mod2000,2004) %>%
  update_model_id()
edit_model(mod2004)

mod2005 <- copy_model_from(mod2000,2005) %>%
  update_model_id()
edit_model(mod2005)

submit_model(
  mod2001,
  .bbi_args = list(overwrite = TRUE, threads = 2)
)

submit_model(
  mod2002,
  .bbi_args = list(overwrite = TRUE, threads = 2)
)

submit_model(
  mod2003,
  .bbi_args = list(overwrite = TRUE, threads = 2)
)

submit_model(
  mod2004,
  .bbi_args = list(overwrite = TRUE, threads = 2)
)

submit_model(
  mod2005,
  .bbi_args = list(overwrite = TRUE, threads = 2)
)

