### Purpose -------------------------------
##' Produce diagnostic plots for the report using an Rmd template 
##' (located in the diagnostics-templates-pk folder)
##' By defining model specific characteristics here, we can simply render the 
##' Rmd template to produce the diagnostics relevant to a specific model
##'
##' Support for parameterized reports can be found 
##' https://bookdown.org/yihui/rmarkdown/parameterized-reports.html
##' 

### Libraries ----------------------------
library(tidyverse)
library(glue)
library(bbr)
library(here)
library(RhpcBLASctl)

### Script ----------------------------
thisScript <- "pk-model-diagnostics-report.R"

### Model directory -------------------
modelDir <- here("model/pk")

# This script includes a `run_sims()` function for generating simulation-based
# diagnostics, and writing them to an RDS file:
#   - EPRED (median/mean and percentiles)
#   - IPRED (median/mean and percentiles, if .iph files available)
#   - NPDE (using median/mean of posterior, or full posterior)
#   - EWRES (using median/mean of posterior, or full posterior)
# It also replaces NONMEM-output ETAs with medians of posterior ETAs from .iph
# files if available.
#
# These simulations require a path to an mrgsolve model.
#
# Although these simulation-based diagnostics are available from NONMEM table
# files for each individual chain, for any final model run with BAYES/NUTS it is
# recommended that these simulations are run using the full posterior using the
# methods implemented in this code.
source(here("script/run-sims-npde.R"))

source(here("script/functions-diagnostics-npde.R"))

# Parallel options --------------------------------------------------------
# Ensure that matrix algebra is not threaded when parallelising simulations 
omp_set_num_threads(1)
blas_set_num_threads(1)
blas_get_num_procs()

# Check number of cores
options(future.fork.enable = TRUE)
future::plan(future::multicore, workers = parallelly::availableCores() - 1)
opt <- furrr::furrr_options(seed = TRUE)
  

# Please go to the
# template-bayes-report.Rmd to see the fields that can be edited directly in this R
# script.

### run 1000: Demo example----------------------------
modelName <- "1000"


### Run this code to regenerate the simulations
system.time({
  run_sims_npde(
    mod_bbr = read_model(file.path(modelDir, modelName, glue("{modelName}_1"))),
    mrgsolve_path = here(glue("model/mrgsolve/{modelName}.mod")),
    out_path = file.path(modelDir, modelName, glue("diag-sims-{modelName}.rds")),
    log_dv = FALSE,
    n_post = 2000
  )
})

### Run this code to regenerate the diagnostics
system.time({
  rmarkdown::render(
    here("script/diagnostic-templates/template-bayes-report.Rmd"),
    params = list(
      run = modelName,
      logDV = FALSE,
      # Specify the path of the simulation output file
      sims_output_path = file.path(modelDir, modelName,
                                   glue("diag-sims-{modelName}.rds"))
    ),
    output_dir = file.path(modelDir, modelName),
    output_file = glue("diagnostic-plots-{modelName}.html")
  )
})

utils::browseURL(file.path(modelDir, modelName,
                           glue("diagnostic-plots-{modelName}.html")))


####Pediatric example#####

modelName <- "2000"

# After an initial look at the diagnostics, run the `run-sims()`
# function to generate the appropriate report diagnostics.

### Run this code to regenerate the simulations - only needed if model re-run
system.time({
  run_sims_npde(
    mod_bbr = read_model(file.path(modelDir, modelName, glue("{modelName}_1"))),
    mrgsolve_path = here(glue("model/mrgsolve/{modelName}.mod")),
    out_path = file.path(modelDir, modelName, glue("diag-sims-{modelName}.rds")),
    log_dv = FALSE,
    n_post = 1500
  )
})

### Run this code to regenerate the diagnostics
system.time({
  rmarkdown::render(
    here("script/diagnostic-templates/template-bayes-report-ped.Rmd"),
    params = list(
      run = modelName,
      logDV = FALSE,
      yspec = "atorvWrkShop3.yml",
      n_thin = 1,
      contCov = c('WT','NUM'),
      catCov = c("FORM","DOSE"),
      etas = c("ETA1//ETA-CL", "ETA2//ETA-V2/F"),
      # Specify the path of the simulation output file
      sims_output_path = file.path(modelDir, modelName,
                                   glue("diag-sims-{modelName}.rds"))
    ),
    output_dir = file.path(modelDir, modelName),
    output_file = glue("diagnostic-plots-{modelName}.html")
  )
})

utils::browseURL(file.path(modelDir, modelName,
                           glue("diagnostic-plots-{modelName}.html")))




