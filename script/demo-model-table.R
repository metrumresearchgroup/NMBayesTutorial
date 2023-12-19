#' ---
#' title: Param tables
#' ---
#'


# load package ------------------------------------------------------------

library(tidyverse)
library(pmtables)
library(bbr)
library(bbr.bayes)
library(yspec)
library(glue)
library(here)
library(posterior)

# directory ---------------------------------------------------------------

runno <- "1000"

params <- list(
  run = runno,
  script = "demo-model-table.R",
  # Directory for the model output
  modelDir = here("model/pk"),
  # Output directory for the derived tables
  tables = here(glue("deliv/table/pk/{runno}")),
  # Analysis dataset
  data = here("data/derived/analysis3.csv"),
  # Key with parameter metadata
  key = here(glue("script/parameter-key-{runno}.yaml"))
)

if (!file.exists(params$tables)) dir.create(params$tables)


# proj and data dir set up
projectDir <- here()
dataDir <- file.path(projectDir, "data")
sourceDataDir <- file.path(dataDir, "source")
derivedDataDir <- file.path(dataDir, "derived")

thisScript <- "demo-model-table.R"


# helper functions --------------------------------------------------------

source(here("script/functions-table.R"))


# read outputs ------------------------------------------------------------

mod_bbr <- read_model(file.path(params$modelDir, params$run))
draws <- read_fit_model(mod_bbr)

draws_param <- draws %>% 
  subset_draws(variable = c("THETA", "OMEGA", "SIGMA"))

# Calculate shrinkage from post hoc ETAs, or from the .shk files if the .iph
# files do not exist
shk0 <- shrinkage(mod_bbr)
omegas <- variables(draws) %>% 
  str_subset("^OMEGA") %>% 
  # diagonals only
  str_subset("\\[(\\d+),\\1\\]")
shk <- tibble(parameter = omegas, shrinkage = shk0)

# parameter estimates -----------------------------------------------------

ptable <- summarize_draws(
  draws_param,
  mean,
  median,
  ~ quantile2(.x, probs = c(0.025, 0.975)),
  rhat,
  ess_bulk,
  ess_tail
) %>% 
  filter(!is.na(rhat)) %>%  # only include non-fixed parameters
  rename(
    name = variable,
    qlo = q2.5,
    qhi = q97.5
  )

param_key <- yaml_as_df(params$key) %>% 
  rename(name = .row)

param_df <- param_key %>%
  left_join(ptable) %>%
  mutate(across(c(mean, median, qlo, qhi), function(.x) {
    case_when(
      trans == "logTrans" ~ exp(.x),
      trans == "lognormalOm" ~ sqrt(exp(.x) - 1) * 100,
      trans == "propErr" ~ sqrt(.x) * 100,
      trans == "addErr" ~ sqrt(.x),
      TRUE ~ .x
    )
  })) %>%
  mutate(
    num = str_extract(name, "[0-9,]+"),
    greek = case_when(
      trans == "logTrans" & panel %in% c("struct", "cov") ~ glue("\\exp(\\theta_{<<num>>})", .open = "<<", .close = ">>"),
      # trans == "invlogit_1" & panel %in% c("struct", "cov") ~ glue("1 / (1 + \\exp(-\\theta_{<<num>>})) \\times 100\\%", .open = "<<", .close = ">>"),
      trans == "none" & panel %in% c("struct", "cov") ~ glue("\\theta_{<<num>>}", .open = "<<", .close = ">>"),
      trans == "lognormalOm" & panel %in% c("iiv") ~ glue("\\Omega_{<<num>>}", .open = "<<", .close = ">>"),
      trans == "none" & panel %in% c("iiv") ~ glue("\\Omega_{<<num>>}", .open = "<<", .close = ">>"),
      trans %in% c("propErr") ~ glue("\\Sigma_{<<num>>}", .open = "<<", .close = ">>"),
      trans %in% c("addErr") ~ glue("\\Sigma_{<<num>>}", .open = "<<", .close = ">>"),
      TRUE ~ NA_character_
    ),
    greek = mathMode(greek)
  ) %>%
  left_join(shk %>% rename(name = parameter, shk = shrinkage)) %>%
  mutate(across(c(mean, median, qlo, qhi, rhat, shk), ~ sig(.x, maxex = 4))) %>%
  mutate(across(c(ess_bulk, ess_tail), round)) %>%
  mutate(across(where(is.character), function(.x) {
    if_else(str_trim(.x) == "NA", glue(""), .x)
  })) %>%
  mutate(ci = glue("({qlo}, {qhi})")) %>%
  select(panel, abb, greek, desc, median, ci, ess_bulk, ess_tail, rhat, shk)

print(param_df, n = Inf)

# Define footnotes --------------------------------------------------------

footerSource <- paste0("source: ", thisScript)
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
  mutate(panel = case_when(
    panel == "struct" ~ glue("Structural model"),
    panel == "cov" ~ glue("Covariate effects"),
    panel == "iiv" ~ glue("Interindividual variability"),
    panel == "propErr" ~ glue("Residual variability"),
    TRUE ~ panel
  )) %>%
  stable(
    panel = "panel",
    align = cols_center(
      desc = col_ragged(7),
      abb = "l"
    ),
    cols_blank = c("abb", "greek", "desc"),
    cols_rename = c(
      "Median" = "median",
      "95\\% CDI" = "ci",
      "Bulk ESS" = "ess_bulk",
      "Tail ESS" = "ess_tail",
      "$\\hat{R}$" = "rhat",
      "Shrinkage (\\%)" = "shk"
    ),
    r_file = thisScript,
    output_file = file.path(params$tables, glue("{run}-param-tab.tex")),
    note_config = noteconf(type = "minipage", width = 1),
    notes = c(footLog, footAbbrev, footDerive1, footDerive2, footDerive3, footDerive4)
  ) %>%
  as_lscape()

# struct_pm %>% st2report() #Check output
# struct_pm %>% stable_save() #Save output

# Table pdf rendering
st2doc(
  struct_pm,
  output_dir = params$tables,
  output_file = glue("{run}_preview_param_table_1.pdf"),
  landscape = TRUE
)


#### Make parameter table for pediatric PK example####

runno <- "2000"

params <- list(
  run = runno,
  n_chain = 4,
  script = "demo-model-table.R",
  modelDir = here("model/pk"),
  tables = here(glue("deliv/table/pk/{runno}")),
  data = here("data/derived/atorvWrkShop2.csv"),
  key = here(glue("script/parameter-key-{runno}.yaml"))
)

if (!file.exists(params$tables)) dir.create(params$tables)


# read outputs ------------------------------------------------------------

mod_bbr <- read_model(file.path(params$modelDir, params$run))
draws <- read_fit_model(mod_bbr)

draws_param <- draws %>% 
  subset_draws(variable = c("THETA", "OMEGA", "SIGMA"))

# Calculate shrinkage from post hoc ETAs, or from the .shk files if the .iph
# files do not exist
shk0 <- shrinkage(mod_bbr)
omegas <- variables(draws) %>% 
  str_subset("^OMEGA") %>% 
  # diagonals only
  str_subset("\\[(\\d+),\\1\\]")
shk <- tibble(parameter = omegas, shrinkage = shk0)

# parameter estimates -----------------------------------------------------

ptable <- summarize_draws(
  draws_param,
  mean,
  median,
  ~ quantile2(.x, probs = c(0.025, 0.975)),
  rhat,
  ess_bulk,
  ess_tail
) %>% 
  filter(!is.na(rhat)) %>%  # only include non-fixed parameters
  rename(
    name = variable,
    qlo = q2.5,
    qhi = q97.5
  )

param_key <- yaml_as_df(params$key) %>% 
  rename(name = .row)

# param_df

param_df <- param_key %>%
  left_join(ptable) %>%
  mutate(across(c(mean, median, qlo, qhi), function(.x) {
    case_when(
      trans == "logTrans" ~ exp(.x),
      trans == "lognormalOm" ~ sqrt(exp(.x) - 1) * 100,
      trans == "propErr" ~ sqrt(.x) * 100,
      trans == "addErr" ~ sqrt(.x),
      TRUE ~ .x
    )
  })) %>%
  mutate(
    num = str_extract(name, "[0-9,]+"),
    greek = case_when(
      trans == "logTrans" & panel %in% c("struct", "cov") ~ glue("\\exp(\\theta_{<<num>>})", .open = "<<", .close = ">>"),
      # trans == "invlogit_1" & panel %in% c("struct", "cov") ~ glue("1 / (1 + \\exp(-\\theta_{<<num>>})) \\times 100\\%", .open = "<<", .close = ">>"),
      trans == "none" & panel %in% c("struct", "cov") ~ glue("\\theta_{<<num>>}", .open = "<<", .close = ">>"),
      trans == "lognormalOm" & panel %in% c("iiv") ~ glue("\\Omega_{<<num>>}", .open = "<<", .close = ">>"),
      trans == "none" & panel %in% c("iiv") ~ glue("\\Omega_{<<num>>}", .open = "<<", .close = ">>"),
      trans %in% c("propErr") ~ glue("\\Sigma_{<<num>>}", .open = "<<", .close = ">>"),
      trans %in% c("addErr") ~ glue("\\Sigma_{<<num>>}", .open = "<<", .close = ">>"),
      TRUE ~ NA_character_
    ),
    greek = mathMode(greek)
  ) %>%
  left_join(shk %>% rename(name = parameter, shk = shrinkage)) %>%
  mutate(across(c(mean, median, qlo, qhi, rhat, shk), ~ sig(.x, maxex = 4))) %>%
  mutate(across(c(ess_bulk, ess_tail), round)) %>%
  mutate(across(where(is.character), function(.x) {
    if_else(str_trim(.x) == "NA", glue(""), .x)
  })) %>%
  mutate(ci = glue("({qlo}, {qhi})")) %>%
  select(panel, abb, greek, desc, median, ci, ess_bulk, ess_tail, rhat, shk)


print(param_df, n = Inf)


# Define footnotes --------------------------------------------------------

footerSource <- paste0("source: ", thisScript)
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
  mutate(panel = case_when(
    panel == "struct" ~ glue("Structural model"),
    panel == "cov" ~ glue("Covariate effects"),
    panel == "iiv" ~ glue("Interindividual variability"),
    panel == "propErr" ~ glue("Residual variability"),
    TRUE ~ panel
  )) %>%
  stable(
    panel = "panel",
    align = cols_center(
      desc = col_ragged(7),
      abb = "l"
    ),
    cols_blank = c("abb", "greek", "desc"),
    cols_rename = c(
      "Median" = "median",
      "95\\% CDI" = "ci",
      "Bulk ESS" = "ess_bulk",
      "Tail ESS" = "ess_tail",
      "$\\hat{R}$" = "rhat",
      "Shrinkage (\\%)" = "shk"
    ),
    r_file = thisScript,
    output_file = file.path(params$tables, glue("{run}-param-tab.tex")),
    note_config = noteconf(type = "minipage", width = 1),
    notes = c(footLog, footAbbrev, footDerive1, footDerive2, footDerive3, footDerive4)
  ) %>%
  as_lscape()

# struct_pm %>% st2report()
# struct_pm %>% stable_save()

# Table pdf
st2doc(
  struct_pm,
  output_dir = params$tables,
  output_file = glue("TableS1.pdf"),
  landscape = TRUE
)
