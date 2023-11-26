### Applied simulations using full posterior distribution###

#' # Required packages
library(tidyverse)
library(glue)
library(bbr)
library(bbr.bayes)
library(here)
library(withr)
library(RhpcBLASctl)
library(yspec)
library(mrgsolve)
library(mrggsave)
library(vpc)
library(future)
library(future.apply)
library(purrr)
library(furrr)
library(mrgmisc)
library(posterior)

options(mrggsave.dir = here("deliv/figure"), mrg.script = "applied-sims.R")
theme_set(pmplots::pm_theme())
plan(multisession)
set.seed(01262022)
opt <- furrr_options(seed = TRUE)

# source(here("script/functions-diagnostics-npde.R"))

# Set up directories and model numbers to pull in posterior
out_data_dir <- here("script", "simout")
modelDir <- here("model/pk")
modelName <- "2000"
mod_bbr <- read_model(file.path(modelDir, modelName))
mrgsolve_path <- here(glue("model/mrgsolve/{modelName}.mod"))

n_post <- 1000

draws <- read_fit_model(mod_bbr)
n_chains <- nchains(draws)
n_iter <- niterations(draws)
n_draws <- ndraws(draws)

draws_param <- draws %>% 
  subset_draws(variable = c("THETA", "OMEGA", "SIGMA")) %>% 
  as_draws_df() %>% 
  # e.g., THETA[1] -> THETA1
  rename_with(~ str_remove_all(.x, "\\[|\\]"))
draws_param

withr::with_seed(01242022, {
  post <- draws_param %>%
    slice_sample(n = n_post) %>%
    mutate(ITERATION = row_number())
})


# Read in mrgsolve model
mod_mrgsolve <- mread(mrgsolve_path)
param(mod_mrgsolve)
omat(mod_mrgsolve)
smat(mod_mrgsolve)

# Set up "population" dosing dataset
# Use analysis dataset to generate a template population

adset <- read_csv(file = here("data/derived/atorvWrkShop3.csv"), na = ".")

withr::with_seed(11012022, {
  popdat <- adset %>%
    distinct(ID, WT) %>%
    resample_df(
      key_cols = "ID",
      n = 1000
    ) %>%
    select(ID, WT) %>%
    mutate(ID = 1:1000)
})


RUN_SIMS <- FALSE
# Set up iterative simulation
if (isTRUE(RUN_SIMS)) {
  withr::with_seed(01242022, {
    # NOTE: I'm just using arbitrary doses here for means of demonstration
    map(c(1000, 750, 500, 250), function(xamt) {
      simres <- future_map_dfr(
        chunk_df(post, ITERATION, .nchunks = 16),
        .options = opt,
        function(.chunk) {
          .chunk <- ungroup(.chunk)
          map_dfr(
            .chunk$ITERATION,
            function(iter, params = .chunk, simmod = mod_mrgsolve, Xamt = xamt,
                     Popdat = popdat) {
              uc <- filter(params, ITERATION == iter)

              ddat <- popdat %>%
                mutate(
                  EVID = 1,
                  TIME = 0,
                  CMT = 1,
                  AMT = (Xamt / 558.64) * 1000,
                  SS = 1,
                  ADDL = 2,
                  II = 24,
                  DOSE = Xamt
                )

              # Set up tgrid
              simt <- tgrid(
                start = 0,
                end = 48,
                delta = 48
              )


              out <- simmod %>%
                data_set(ddat) %>%
                param(uc %>% select(starts_with("THETA"))) %>%
                omat(uc %>% as_bmat("OMEGA")) %>%
                smat(uc %>% as_dmat("SIGMA")) %>%
                carry_out(DOSE) %>%
                mrgsim_df(obsonly = TRUE, tgrid = simt) %>%
                filter(TIME == 48) %>%
                select(ID, Y, DOSE) %>%
                mutate(run = iter) %>%
                group_by(DOSE, run) %>%
                summarize(perc = sum(Y < 60) / n() * 100, .groups = "drop")


              return(out)
            }
          )
        }
      )

      saveRDS(simres, file = here(
        out_data_dir,
        paste("simresPop", xamt, ".RDS",
          sep = ""
        )
      ))
      rm(simres)
      gc()
    })
  })
}

# Read in simulated output.
simres2 <- readRDS(file = here(
  out_data_dir,
  "simresPop250.RDS"
))
simres3 <- readRDS(file = here(
  out_data_dir,
  "simresPop500.RDS"
))
simres4 <- readRDS(file = here(
  out_data_dir,
  "simresPop750.RDS"
))
simres5 <- readRDS(file = here(
  out_data_dir,
  "simresPop1000.RDS"
))

# Stack the datasets
simres <- bind_rows(
  simres2, simres3,
  simres4, simres5
) %>%
  mutate(POP = case_when(
    DOSE == 250 ~ "Dose 1",
    DOSE == 500 ~ "Dose 2",
    DOSE == 750 ~ "Dose 3",
    DOSE == 1000 ~ "Dose 4",
    TRUE ~ "CHECK"
  ))

simressum <- simres %>%
  group_by(POP) %>%
  summarize(qlo = quantile(perc, 0.1), .groups = "drop")
simressum

# Generate plot
p1 <- ggplot(simres, aes(x = perc, group = POP)) +
  geom_density(aes(fill = POP), alpha = 0.55) +
  scale_color_brewer(palette = "Set1") +
  labs(
    y = "Density",
    x = "Percent of Population Below Cmin Threshold",
    fill = "Dose", color = "Dose",
    caption = "NOTE: Distributions of percentages are derived from 1000 posterior simulations of 1000 subjects.
       The dashed horizotal line is the threshold over which 90% of the simulations are desired to exceed."
  ) +
  theme(legend.position = "top") +
  guides(fill = guide_legend(ncol = 3)) +
  geom_vline(xintercept = 90, lty = 2, color = "steelblue")
p1

mrggsave(p1, tag = "appliedsim", height = 7.5, width = 7.5)
