#' --- 
#' title: Param tables
#' ---
#' 


# load package ------------------------------------------------------------

library(tidyverse)
library(pmtables)
library(bbr)
library(yspec)
library(glue)
library(here)

# directory ---------------------------------------------------------------

setwd(here())

runno <- '1000'

params <- list(
  run = runno,
  #number of chains that were run
  n_chain = 4,
  script = "demo-model-table.R",
  #Directory for the model output
  modelDir = here::here("model/pk"),
  #Output directory for the derived tables
  tables = here::here(glue("deliv/table/pk/{runno}")), 
  #Analsis dataset
  data = "../data/derived/analysis3.csv"
)

if (!file.exists(params$tables)) dir.create(params$tables)


# proj and data dir set up
projectDir <- here::here()
dataDir <- file.path(projectDir, "data")
sourceDataDir <- file.path(dataDir, "source")
derivedDataDir <- file.path(dataDir, "derived")

thisScript <- "demo-model-table.R"


# helper functions --------------------------------------------------------

source(here::here("script/functions-diagnostics.R"))
source(here::here("script/functions-table.R"))


# read outputs ------------------------------------------------------------

# Read in parameter values from ext files (posterior samples)
ext0 <- map_dfr(seq_len(params$n_chain), function(.chain) {
  data.table::fread(
    file = file.path(params$modelDir, params$run, glue("{params$run}_{.chain}"),
                     glue("{params$run}_{.chain}.ext"))
  ) %>% 
    as_tibble() %>% 
    mutate(chain = .chain)
})

# remove early iterations if necessary
ext <- ext0 %>%  filter(ITERATION > 0)

# ext table
ext_tbl <- ext %>% select(-ITERATION, -MCMCOBJ)
n_param <- ncol(ext_tbl) - 1  # don't include chain
n_iter <- nrow(ext) / params$n_chain

#Generate param_array object from NONMEM output for use with Rstan functions
param_array <- array(
  double(n_iter * params$n_chain * n_param),
  dim = c(n_iter, params$n_chain, n_param),
  dimnames = list(NULL, NULL, setdiff(names(ext_tbl), "chain"))
)
for (.chain in seq_len(params$n_chain)) {
  param_array[,.chain,] <- ext_tbl %>% 
    filter(chain == .chain) %>% 
    select(-chain) %>% 
    as.matrix()
}

##summarize posteriors
ext_med <- ext %>% summarize_all(median)
ext_mean <- ext %>% summarize_all(mean)
ext_lo <- ext %>% summarize_all(quantile, prob = 0.025)
ext_hi <- ext %>% summarize_all(quantile, prob = 0.975)

# Read iph files (individual posteriors samples for ETAs)
fn_iph <- file.path(
  params$modelDir,
  params$run,
  glue("{params$run}_{params$n_chain}"),
  glue("{params$run}_{params$n_chain}.iph")
)

if (file.exists(fn_iph)) {
  
  iph <- get_chain_files(params$modelDir, params$run, params$n_chain, "iph") %>% 
    filter(ITERATION > 0)
  ipar <- iph %>% select(chain, ITERATION, ID, starts_with("ETA"))
  iobj <- iph %>% select(chain, ITERATION, ID, MCMCOBJ)
  
  # calculate shrinkage
  shk0 <- baysh(ipar, ext)
  
} else {
  
  shk0 <- get_chain_files(params$modelDir, params$run, params$n_chain, "shk") %>% 
    # From Intro to NM 7: "Type 4=%Eta shrinkage SD version"
    filter(TYPE == 4) %>% 
    select(starts_with("ETA")) %>% 
    summarise(across(.fns = median)) %>% 
    pivot_longer(everything(), names_to = "eta", values_to = "shrinkage")
  
}

shk <- shk0 %>% 
  mutate(
    shrinkage = sig(shrinkage),
    eta = str_extract(eta, "[\\d]+"),
    parameter = glue("OMEGA({eta},{eta})")
  ) %>% 
  select(-eta)

# parameter estimates -----------------------------------------------------

# 1st chain details 
mod <- read_model(file.path(params$modelDir, params$run, glue("{params$run}_1")))
sum <- mod %>% model_summary()

n_theta <- length(sum$parameter_names$theta)
n_omega <- length(sum$parameter_names$omega)
n_sigma <- length(sum$parameter_names$sigma)
size_omega <- (-1 + sqrt(1 + 4*1*n_omega*2)) / 2
size_sigma <- (-1 + sqrt(1 + 4*1*n_sigma*2)) / 2
fix_diag <- tibble(
  parameter = sum$parameter_names %>% do.call(c, .),
  fixed = sum$parameters_data[[1]]$fixed %>% do.call(c, .) %>% as.logical,
  diag = c(rep(NA, n_theta), block(size_omega), block(size_sigma))
)

# parameter table generation using rstan functions

ptable <- rstan::monitor(param_array, warmup = 0, print = FALSE) %>% 
  as.matrix %>%
  as.data.frame() %>% 
  mutate(
    parameter = rownames(.)
    #mean_med_diff = (mean - `50%`) / `50%` * 100
  ) %>% 
  left_join(fix_diag) %>% 
  filter(!fixed)

ptable %>%
  rename(pct2.5 = "2.5%", median = "50%", pct97.5 = "97.5%", Neff = "n_eff") %>%
  mutate("95% CI" = glue("({pct2.5},{pct97.5})")) %>% 
  select(parameter, mean, median, "95% CI", Bulk_ESS, Tail_ESS, Rhat) %>%
  left_join(shk, by = "parameter") %>% 
  knitr::kable(caption = "Summary of Model Parameter Estimates.") 

# parameter key labels with Latex formatting
param_key <- tribble(
  ~name, ~abbr, ~desc, ~trans, ~type,
  "THETA1",     "$\\text{KA}$ (1/hr)",          "First order absorption rate constant",       "exp",         "struct", 
  "THETA2",     "$\\text{V2/F}$ (L)",           "Apparent central volume of distribution",    "exp",         "struct",
  "THETA3",     "$\\text{CL/F}$ (L/hr)",        "Apparent clearance",                         "exp",         "struct",
  "THETA4",     "$\\text{V3/F}$ (L)",           "Apparent peripheral volume of distribution", "exp",         "struct", 
  "THETA5",     "$\\text{Q/F}$ (L/hr)",         "Apparent intercompartmental clearance",      "exp",         "struct",
  "OMEGA(1,1)",    "$\\Omega_\\mathrm{KA}$",       "IIV-KA (CV\\%)",                           "lognormalOm", "iiv",
  "OMEGA(2,2)",    "$\\Omega_\\mathrm{V2/F}$",       "IIV-V2/F (CV\\%)",                           "lognormalOm", "iiv",
  "OMEGA(3,3)",    "$\\Omega_\\mathrm{CL/F}$",       "IIV-CL/F (CV\\%)",                           "lognormalOm", "iiv",
  "OMEGA(2,1)", "V2/F-KA", "Covariance of V2/F - KA", "none", "iiv",
  "OMEGA(3,1)", "CL/F-KA", "Covariance of CL/F - KA", "none", "iiv",
  "OMEGA(3,2)", "CL/F-V2/F", "Covariance of CL/F - V2/F", "none", "iiv",
  "SIGMA(1,1)", "$\\Sigma_\\mathrm{11}$",       "RUV - proportional (CV\\%)",   "propErr",     "RV"
)


# param_df formatting for Latex aesthetics. 

param_df <- param_key %>% 
  left_join(
    ptable %>% 
      as_tibble() %>% 
      select(
        name = parameter, mean, median = Q50, qlo = `2.5%`, qhi = `97.5%`,
        Bulk_ESS, Tail_ESS, Rhat
      )
  ) %>% 
  mutate(across(c(mean, median, qlo, qhi), function(.x) {
    case_when(
      trans == "exp" ~ exp(.x),
      trans == "lognormalOm" ~ sqrt(exp(.x)-1)*100,
      trans == "propErr" ~ sqrt(.x)*100,
      trans == "addErr" ~ sqrt(.x),
      TRUE ~ .x
    )
  })) %>% 
  mutate(
    num = str_extract(name, "[0-9,]+"),
    greek = case_when(
      trans == "exp" & type %in% c("struct", "cov") ~ glue("\\exp(\\theta_{<<num>>})", .open = "<<", .close = ">>"),
      # trans == "invlogit_1" & type %in% c("struct", "cov") ~ glue("1 / (1 + \\exp(-\\theta_{<<num>>})) \\times 100\\%", .open = "<<", .close = ">>"),
      trans == "none" & type %in% c("struct", "cov") ~ glue("\\theta_{<<num>>}", .open = "<<", .close = ">>"),
      trans == "lognormalOm" & type %in% c("iiv") ~ glue("\\Omega_{<<num>>}", .open = "<<", .close = ">>"),
      trans == "none" & type %in% c("iiv") ~ glue("\\Omega_{<<num>>}", .open = "<<", .close = ">>"),
      trans %in% c("propErr") ~ glue("\\Sigma_{<<num>>}", .open = "<<", .close = ">>"),
      trans %in% c("addErr") ~ glue("\\Sigma_{<<num>>}", .open = "<<", .close = ">>"),
      TRUE ~ NA_character_
    ),
    greek = mathMode(greek)
  ) %>% 
  left_join(shk %>% rename(name = parameter, shk = shrinkage)) %>% 
  mutate(across(c(mean, median, qlo, qhi, Rhat), ~sig(.x,maxex=4))) %>% 
  mutate(across(where(is.character), function(.x) {
    if_else(str_trim(.x) == "NA", glue(""), .x)
  })) %>% 
  mutate(ci = glue("({qlo}, {qhi})")) %>% 
  select(type, abbr, greek, desc, median, ci, Bulk_ESS, Tail_ESS, Rhat, shk)

print(param_df, n = Inf)


# Define footnotes --------------------------------------------------------

footerSource <- paste0('source: ', thisScript)
footAbbrev <- "Abbreviations: CDI = credible interval;
  ESS = effective sample size;
  $\\hat{R}$ = Gelman-Rubin diagnostic;
  IIV = inter-individual variability;
  RV = residual variability;
  CV = coefficient of variation; 
  SD = standard deviation"
footAbbrev2 <- "Abbreviations: CDI = credible interval;
CI = confidence interval;
IIV = inter-individual variability;
IOV = inter-occasion variability;
RV = residual variability"
footAbbrev_Om <- "Abbreviations: CDI = credible interval;
  ESS = effective sample size;
  $\\hat{R}$ = Gelman-Rubin diagnostic"
footDerive1 <- "Credible intervals calculated from Bayesian posteriors"
footDerive2 <- "CV\\% of omegas = sqrt(exp(estimate) - 1) * 100"
footDerive3 <- "CV\\% of sigma = sqrt(estimate) * 100"
footDerive4 <- "SD of sigma = sqrt(estimate)"
footLog <- "Parameters estimated in the log-domain were back-transformed for clarity"

run <- params$run


# Create table ------------------------------------------------------------

##  create table with structural parameters and RV using pmtables package
struct_pm <- param_df %>% 
  mutate(type = case_when(
    type == "struct" ~ glue("Structural model"),
    type == "cov" ~ glue("Covariate effects"),
    type == "iiv" ~ glue("Interindividual variability"),
    type == "propErr" ~ glue("Residual variability"),
    TRUE ~ type
  )) %>%
  stable(
    panel = "type",
    align = cols_center(desc = col_ragged(7),
                        abbr = "l"),
    cols_blank = c("abbr", "greek", "desc"),
    cols_rename = c(
      "Median" = "median", 
      "95\\% CDI" = "ci",
      "Bulk ESS" = "Bulk_ESS",
      "Tail ESS" = "Tail_ESS",
      "$\\hat{R}$" = "Rhat",
      "Shrinkage (\\%)" = "shk"
    ),
    r_file = thisScript,
    output_file = file.path(params$tables, glue("{run}-param-tab.tex")),
    note_config = noteconf(type = "minipage", width = 1),
    notes = c(footLog, footAbbrev, footDerive1, footDerive2, footDerive3, footDerive4)
  ) %>% as_lscape()

# struct_pm %>% st2report() #Check output
# struct_pm %>% stable_save() #Save output

# Table pdf rendering
st2doc(
  struct_pm,
  output_dir = params$tables,
  output_file = glue("{run}_preview_param_table_1.pdf"),
  landscape = TRUE
)


####Make parameter table for pediatric PK example####

runno <- '2000'

params <- list(
  run = runno, 
  n_chain = 4,
  script = "demo-model-table.R",
  modelDir = here::here("model/pk"),
  tables = here::here(glue("deliv/table/pk/{runno}")), 
  data = "../data/derived/atorvWrkShop2.csv"
)

if (!file.exists(params$tables)) dir.create(params$tables)


# read outputs ------------------------------------------------------------

# Read in parameter values from ext files 
ext0 <- map_dfr(seq_len(params$n_chain), function(.chain) {
  data.table::fread(
    file = file.path(params$modelDir, params$run, glue("{params$run}_{.chain}"),
                     glue("{params$run}_{.chain}.ext"))
  ) %>% 
    as_tibble() %>% 
    mutate(chain = .chain)
})

# remove early iterations if necessary
ext <- ext0 %>%  filter(ITERATION > 0)

# ext table
ext_tbl <- ext %>% select(-ITERATION, -MCMCOBJ)
n_param <- ncol(ext_tbl) - 1  # don't include chain
n_iter <- nrow(ext) / params$n_chain
param_array <- array(
  double(n_iter * params$n_chain * n_param),
  dim = c(n_iter, params$n_chain, n_param),
  dimnames = list(NULL, NULL, setdiff(names(ext_tbl), "chain"))
)
for (.chain in seq_len(params$n_chain)) {
  param_array[,.chain,] <- ext_tbl %>% 
    filter(chain == .chain) %>% 
    select(-chain) %>% 
    as.matrix()
}

##summarize posteriors
ext_med <- ext %>% summarize_all(median)
ext_mean <- ext %>% summarize_all(mean)
ext_lo <- ext %>% summarize_all(quantile, prob = 0.025)
ext_hi <- ext %>% summarize_all(quantile, prob = 0.975)

# Read iph files 
fn_iph <- file.path(
  params$modelDir,
  params$run,
  glue("{params$run}_{params$n_chain}"),
  glue("{params$run}_{params$n_chain}.iph")
)

if (file.exists(fn_iph)) {
  
  iph <- get_chain_files(params$modelDir, params$run, params$n_chain, "iph") %>% 
    filter(ITERATION > 0)
  ipar <- iph %>% select(chain, ITERATION, ID, starts_with("ETA"))
  iobj <- iph %>% select(chain, ITERATION, ID, MCMCOBJ)
  
  # calculate shrinkage
  shk0 <- baysh(ipar, ext)
  
} else {
  
  shk0 <- get_chain_files(params$modelDir, params$run, params$n_chain, "shk") %>% 
    # From Intro to NM 7: "Type 4=%Eta shrinkage SD version"
    filter(TYPE == 4) %>% 
    select(starts_with("ETA")) %>% 
    summarise(across(.fns = median)) %>% 
    pivot_longer(everything(), names_to = "eta", values_to = "shrinkage")
  
}

shk <- shk0 %>% 
  mutate(
    shrinkage = sig(shrinkage),
    eta = str_extract(eta, "[\\d]+"),
    parameter = glue("OMEGA({eta},{eta})")
  ) %>% 
  select(-eta)

# parameter estimates -----------------------------------------------------

# 1st chain details 
mod <- read_model(file.path(params$modelDir, params$run, glue("{params$run}_1")))
sum <- mod %>% model_summary()

n_theta <- length(sum$parameter_names$theta)
n_omega <- length(sum$parameter_names$omega)
n_sigma <- length(sum$parameter_names$sigma)
size_omega <- (-1 + sqrt(1 + 4*1*n_omega*2)) / 2
size_sigma <- (-1 + sqrt(1 + 4*1*n_sigma*2)) / 2
fix_diag <- tibble(
  parameter = sum$parameter_names %>% do.call(c, .),
  fixed = sum$parameters_data[[1]]$fixed %>% do.call(c, .) %>% as.logical,
  diag = c(rep(NA, n_theta), block(size_omega), block(size_sigma))
)

# parameter table 

ptable <- rstan::monitor(param_array, warmup = 0, print = FALSE) %>% 
  as.matrix %>%
  as.data.frame() %>% 
  mutate(
    parameter = rownames(.)
    #mean_med_diff = (mean - `50%`) / `50%` * 100
  ) %>% 
  left_join(fix_diag) %>% 
  filter(!fixed)

ptable %>%
  rename(pct2.5 = "2.5%", median = "50%", pct97.5 = "97.5%", Neff = "n_eff") %>%
  mutate("95% CI" = glue("({pct2.5},{pct97.5})")) %>% 
  select(parameter, mean, median, "95% CI", Bulk_ESS, Tail_ESS, Rhat) %>%
  left_join(shk, by = "parameter") %>% 
  knitr::kable(caption = "Summary of Model Parameter Estimates.") 

# parameter key
param_key <- tribble(
  ~name, ~abbr, ~desc, ~trans, ~type,
  "THETA1",     "$\\text{CL/F}$ (L/hr)",        "Apparent clearance",                         "exp",         "struct",
  "THETA2",     "$\\text{V2/F}$ (L)",           "Apparent central volume of distribution",    "exp",         "struct",
  "THETA3",     "$\\text{Q/F}$ (L/hr)",         "Apparent intercompartmental clearance",      "exp",         "struct",
  "THETA4",     "$\\text{V3/F}$ (L)",           "Apparent peripheral volume of distribution", "exp",         "struct", 
  "THETA5",     "$\\text{KA}$ (1/hr)",          "First order absorption rate constant",       "exp",         "struct", 
  "OMEGA(1,1)",    "$\\Omega_\\mathrm{CL/F}$",       "IIV-CL/F (CV\\%)",                           "lognormalOm", "iiv",
  "OMEGA(2,1)", "CL/F-V2/F", "Covariance of CL/F - V2/F", "none", "iiv",
  "OMEGA(2,2)",    "$\\Omega_\\mathrm{V2/F}$",       "IIV-V2/F (CV\\%)",                           "lognormalOm", "iiv",
  "SIGMA(1,1)", "$\\Sigma_\\mathrm{11}$",       "RUV - proportional (CV\\%)",   "propErr",     "RV"
)


# param_df 

param_df <- param_key %>% 
  left_join(
    ptable %>% 
      as_tibble() %>% 
      select(
        name = parameter, mean, median = Q50, qlo = `2.5%`, qhi = `97.5%`,
        Bulk_ESS, Tail_ESS, Rhat
      )
  ) %>% 
  mutate(across(c(mean, median, qlo, qhi), function(.x) {
    case_when(
      trans == "exp" ~ exp(.x),
      trans == "lognormalOm" ~ sqrt(exp(.x)-1)*100,
      trans == "propErr" ~ sqrt(.x)*100,
      trans == "addErr" ~ sqrt(.x),
      TRUE ~ .x
    )
  })) %>% 
  mutate(
    num = str_extract(name, "[0-9,]+"),
    greek = case_when(
      trans == "exp" & type %in% c("struct", "cov") ~ glue("\\exp(\\theta_{<<num>>})", .open = "<<", .close = ">>"),
      # trans == "invlogit_1" & type %in% c("struct", "cov") ~ glue("1 / (1 + \\exp(-\\theta_{<<num>>})) \\times 100\\%", .open = "<<", .close = ">>"),
      trans == "none" & type %in% c("struct", "cov") ~ glue("\\theta_{<<num>>}", .open = "<<", .close = ">>"),
      trans == "lognormalOm" & type %in% c("iiv") ~ glue("\\Omega_{<<num>>}", .open = "<<", .close = ">>"),
      trans == "none" & type %in% c("iiv") ~ glue("\\Omega_{<<num>>}", .open = "<<", .close = ">>"),
      trans %in% c("propErr") ~ glue("\\Sigma_{<<num>>}", .open = "<<", .close = ">>"),
      trans %in% c("addErr") ~ glue("\\Sigma_{<<num>>}", .open = "<<", .close = ">>"),
      TRUE ~ NA_character_
    ),
    greek = mathMode(greek)
  ) %>% 
  left_join(shk %>% rename(name = parameter, shk = shrinkage)) %>% 
  mutate(across(c(mean, median, qlo, qhi, Rhat), ~sig(.x,maxex=4))) %>% 
  mutate(across(where(is.character), function(.x) {
    if_else(str_trim(.x) == "NA", glue(""), .x)
  })) %>% 
  mutate(ci = glue("({qlo}, {qhi})")) %>% 
  select(type, abbr, greek, desc, median, ci, Bulk_ESS, Tail_ESS, Rhat, shk)

print(param_df, n = Inf)


# Define footnotes --------------------------------------------------------

footerSource <- paste0('source: ', thisScript)
footAbbrev <- "Abbreviations: CDI = credible interval;
  ESS = effective sample size;
  $\\hat{R}$ = Gelman-Rubin diagnostic;
  IIV = inter-individual variability;
  RV = residual variability;
  CV = coefficient of variation; 
  SD = standard deviation"
footAbbrev2 <- "Abbreviations: CDI = credible interval;
CI = confidence interval;
IIV = inter-individual variability;
IOV = inter-occasion variability;
RV = residual variability"
footAbbrev_Om <- "Abbreviations: CDI = credible interval;
  ESS = effective sample size;
  $\\hat{R}$ = Gelman-Rubin diagnostic"
footDerive1 <- "Credible intervals calculated from Bayesian posteriors"
footDerive2 <- "CV\\% of omegas = sqrt(exp(estimate) - 1) * 100"
footDerive3 <- "CV\\% of sigma = sqrt(estimate) * 100"
footDerive4 <- "SD of sigma = sqrt(estimate)"
footLog <- "Parameters estimated in the log-domain were back-transformed for clarity"

run <- params$run


# Create table ------------------------------------------------------------

##  create table with structural parameters and RV
struct_pm <- param_df %>% 
  mutate(type = case_when(
    type == "struct" ~ glue("Structural model"),
    type == "cov" ~ glue("Covariate effects"),
    type == "iiv" ~ glue("Interindividual variability"),
    type == "propErr" ~ glue("Residual variability"),
    TRUE ~ type
  )) %>%
  stable(
    panel = "type",
    align = cols_center(desc = col_ragged(7),
                        abbr = "l"),
    cols_blank = c("abbr", "greek", "desc"),
    cols_rename = c(
      "Median" = "median", 
      "95\\% CDI" = "ci",
      "Bulk ESS" = "Bulk_ESS",
      "Tail ESS" = "Tail_ESS",
      "$\\hat{R}$" = "Rhat",
      "Shrinkage (\\%)" = "shk"
    ),
    r_file = thisScript,
    output_file = file.path(params$tables, glue("{run}-param-tab.tex")),
    note_config = noteconf(type = "minipage", width = 1),
    notes = c(footLog, footAbbrev, footDerive1, footDerive2, footDerive3, footDerive4)
  ) %>% as_lscape()

# struct_pm %>% st2report()
# struct_pm %>% stable_save()

# Table pdf
st2doc(
  struct_pm,
  output_dir = params$tables,
  output_file = glue("{run}_preview_param_table_ped-pk.pdf"),
  landscape = TRUE
)

