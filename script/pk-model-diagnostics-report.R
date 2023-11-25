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
library(bbr.bayes)
library(here)
library(mrgsolve)

### Script ----------------------------
thisScript <- "pk-model-diagnostics-report.R"

### Model directory -------------------
model_dir <- here("model/pk")
model_dir_sim <- here("model/mrgsolve")

# This script uses `bbr.bayes::nm_join_bayes()` for generating simulation-based
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

#These two scripts contain helper functions that perform simulations and 
#summarizations for the creation of NPDE diagnostics. 
# source(here("script/run-sims-npde.R"))

# source(here("script/functions-diagnostics-npde.R"))

# Parallel options --------------------------------------------------------
# # Ensure that matrix algebra is not threaded when parallelising simulations 
# omp_set_num_threads(1)
# blas_set_num_threads(1)
# blas_get_num_procs()

# Check number of cores
options(future.fork.enable = TRUE)
future::plan(future::multicore, workers = parallelly::availableCores() - 1)
opt <- furrr::furrr_options(seed = TRUE)
  

# Please go to the
# template-bayes-report.Rmd to see the fields that can be edited directly in this R
# script.

### run 1000: Demo example----------------------------
model_name <- "1000"


### Run this code to regenerate the simulations

# Refer to ?nm_join_bayes() for additional details on the function arguments. 

# Decompress .iph files. This is only needed because the file sizes are too
# large for GitHub.
# No guarantees this will work on non-unix OSes
for (i in 1:4) {
  dir_i <- file.path(model_dir, glue("1000/1000-{i}"))
  # zip file was split into multiple parts, so need to join together first
  parts <- list.files(dir_i, glue("1000-{i}.iph.z*"), full.names = TRUE)
  parts <- str_subset(parts, ".iph-all", negate = TRUE)
  new_zip <- file.path(dir_i, glue("1000-{i}.iph-all.zip"))
  system(paste("cat", paste(parts, collapse = " "), ">", new_zip))
  system(paste("unzip", new_zip, "-d", dir_i))
}

set.seed(201)
mod_bbr <- read_model(file.path(model_dir, model_name))
mod_ms <- mread(file.path(model_dir_sim, glue("{model_name}.mod")))
progressr::with_progress({
  sim <- nm_join_bayes(mod_bbr, mod_ms, n_post = 2000)
})
saveRDS(
  sim,
  file.path(model_dir, model_name, glue("diag-sims-{model_name}.rds"))
)

### Run this code to regenerate the diagnostics
system.time({
  rmarkdown::render(
    here("script/diagnostic-templates/template-bayes-report.Rmd"),
    params = list(
      run = model_name,
      n_post = 2000,
      logDV = FALSE,
      # Specify the path of the simulation output file (generated from run_sims_npde)
      sims_output_path = file.path(model_dir, model_name,
                                   glue("diag-sims-{model_name}.rds"))
    ),
    #output directory for the derived diagnostic file
    output_dir = file.path(model_dir, model_name),
    #name of derived diagnostic html
    output_file = glue("diagnostic-plots-{model_name}.html")
  )
})

utils::browseURL(file.path(model_dir, model_name,
                           glue("diagnostic-plots-{model_name}.html")))


####Pediatric example#####

model_name <- "2000"

# After an initial look at the diagnostics, run the `run-sims()`
# function to generate the appropriate report diagnostics.

### Run this code to regenerate the simulations - only needed if model re-run
set.seed(201)
mod_bbr <- read_model(file.path(model_dir, model_name))
mod_ms <- mread(file.path(model_dir_sim, glue("{model_name}.mod")))
progressr::with_progress({
  sim <- nm_join_bayes(mod_bbr, mod_ms, n_post = 1500)
})
saveRDS(
  sim,
  file.path(model_dir, model_name, glue("diag-sims-{model_name}.rds"))
)

### Run this code to regenerate the diagnostics
system.time({
  rmarkdown::render(
    here("script/diagnostic-templates/template-bayes-report-ped.Rmd"),
    params = list(
      run = model_name,
      logDV = FALSE,
      yspec = here("data/spec/atorvWrkShop3.yml"),
      n_thin = 1,
      n_thin2 = 1,
      n_post = 1500,
      contCov = c('WT','NUM'),
      catCov = c("FORM","DOSE"),
      etas = c("ETA1//ETA-CL", "ETA2//ETA-V2/F"),
      # Specify the path of the simulation output file
      sims_output_path = file.path(model_dir, model_name,
                                   glue("diag-sims-{model_name}.rds"))
    ),
    output_dir = file.path(model_dir, model_name),
    output_file = glue("diagnostic-plots-{model_name}.html")
  )
})

utils::browseURL(file.path(model_dir, model_name,
                           glue("diagnostic-plots-{model_name}.html")))




