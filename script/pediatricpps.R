###Prior predictive simulations for pediatric example###

#' # Required packages
library(tidyverse)
library(data.table)
library(glue)
library(bbr)
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


options(mrggsave.dir = here("deliv/figure"), mrg.script = "pediatricpps.R")
theme_set(pmplots::pm_theme())

mrg_vpc_theme = new_vpc_theme(list(
  sim_pi_fill = "steelblue3", sim_pi_alpha = 0.5,
  sim_median_fill = "grey60", sim_median_alpha = 0.5
))


#Read in the mrgsolve model file
modelName <- "2000pps"
modelDir <- here("model/pk")
mrgsolve_path = here(glue("model/mrgsolve/{modelName}.mod"))
mod_mrgsolve <- mread(mrgsolve_path)
param(mod_mrgsolve)
omat(mod_mrgsolve)
smat(mod_mrgsolve)

#Set up "population" dosing dataset
#Use the analysis dataset to generate a template population

adset <- read_csv(file = here('data/derived/atorvWrkShop3.csv'),na='.')

#Check the AMT for a 10 mg dose
adset %>% filter(EVID==1) %>%
  count(AMT, DOSE)

#Set up a hypothetical population dataset to evaluate prior predictive simulation
#1000 subjects, at a dose of 10 mg and we'll look at the PK at SS
#Generate WT by resampling from the distribution in the existing dataset

withr::with_seed(105978665,{
  ddat <- data.frame(ID = 1:1000,
                   EVID = 1,
                   CMT = 1,
                   AMT = 17.9,
                   TIME = 0,
                   II = 24,
                   ADDL = 2,
                   SS = 1,
                   WT = sample(adset %>% distinct(ID,WT) %>% pull(WT),
                               size = 1000,replace = T))
}
)

#Read in parameter vectors to simulated from the priors (2000pps.ctl)
#The control stream was run with $SIML and PRIOR=TRUE
#Note that we need a 1 sample estimation step to get the LKJ priors to be correclty applied

paramvec <- fread(
  file = here('model/pk/2000pps/2000pps.ext'),
  skip=1,
  fill = T) %>% 
  filter(ITERATION==0) %>%
  mutate(samp = 1:n()) %>%
  mutate(across(.fns = as.numeric))


#Simulate across these samples

withr::with_seed(01242022, {
  simres <- map_dfr(1:nrow(paramvec), 
                    function(iter, params = paramvec,
                             simmod = mod_mrgsolve,
                             Popdat = ddat){
    #browser()
    
    uc <- filter(params, samp==iter) 
    
    simt <- tgrid(start = 24,
                  end = 48,
                  delta = 24,
                  add = 24 + c(2,4,6,12))
    
    out <- simmod %>%
      data_set(Popdat) %>%
      param(uc %>% select(starts_with('THETA'))) %>%
      omat(uc %>% as_bmat('OMEGA')) %>%
      smat(uc %>% as_dmat("SIGMA"))  %>%
      mrgsim_df(obsonly=TRUE, tgrid=simt, tad=T) %>%
      mutate(run = iter) %>%
      #Censor BLQ values
      filter(Y>= 0.001)
    
    return(out)
                    }
  )
}
)

#Make a VPC of the simulated data to assess the choices of prior
p1 <- vpc(
  #obs = fsad_mad,
  sim = simres,
  #stratify = "STUDY",
  #obs_cols = list(dv = "DVN"),
  sim_cols=list(dv="Y", sim="run",idv='tad'), 
  log_y = TRUE,
  pi = c(0.05, 0.95),
  ci = c(0.025, 0.975), 
  facet = "columns",
  #show = list(obs_dv = TRUE), 
  vpc_theme = mrg_vpc_theme
)       

#There's a lot of variability, but the interpretation here is that the 
#priors combined with the structural model provide sufficient coverage 
#across the realm of feasibility wrt potential PK concentration profiles
p1
















