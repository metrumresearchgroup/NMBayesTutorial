####Test Comparison Example####
mod3000 <- new_model(file.path(MODEL_DIR, 3000)) 

mod3000 <- read_model(file.path(MODEL_DIR, 3000))

submit_model(
  mod3000,
  .bbi_args = list(overwrite = TRUE),
  .mode = "local"
)

read.table(file.path(MODEL_DIR, "3000/3000.chn"), header = TRUE) %>% 
  select(ITERATION:`OMEGA.3.3.`) %>% 
  pivot_longer(cols = -ITERATION) %>% 
  pivot_wider(names_from = "ITERATION")

run_chains(MODEL_DIR, 3000, .bbi_args = list(
  overwrite = TRUE, parallel = TRUE, threads = 1),
  .mode = "local")

####Compare iph with last posterior sample with table file for ETA values. #####
iph <- data.table::fread(file='../model/pk/3000/3000_1/3000_1.iph', 
                         header = T, skip =1) %>%
  filter(ITERATION > 0 ) %>%
  group_by(ID) %>%
  mutate(across(c(`PHI(1)`:`ETA(5)`),.fns = list(median = median,mean = mean))) %>%
  ungroup() %>%
  distinct(ID, .keep_all = T) %>%
  select(3,17:36)

iphlast <- data.table::fread(file='../model/pk/3000/3000_1/3000_1.iph', 
                         header = T, skip =1) %>%
  filter(ITERATION == max(ITERATION))

phi <- data.table::fread(file='../model/pk/3000/3000_1/3000_1.phi', 
                         header = T, skip =1)

adset <- read_csv(file = here::here('data/derived/atorvWrkShop2.csv'),na='.')
tab <- data.table::fread(file='../model/pk/3000/3000_1/3000.tab', 
                         header = T, skip =1) %>% left_join(adset %>% 
                                                              select(ID,NUM)) %>%
  group_by(ID) %>%
  mutate(across(c(`ETA1`:`ETA5`),.fns = list(median = median,mean = mean))) %>%
  ungroup() %>%
  distinct(ID, .keep_all = T) %>%
  select(ID, ETA1:ETA5) %>%
  arrange(ID)
  


####Proposed structure for supplemental repository######

#Toplevel
#Include a general project repo set up with README
#README: Description on files and instructions on setting up pkgr; link out to "full" Bayes project?
#Make clear all data is simulated

#data
#data/spec: Data specs
#data/derived: .csv format with define documents

#deliv
#deliv/figure: Model specific folders for diagnostics; folder for applied simulation
#deliv/table: Model specific folders for parameter tables

#model
#model/mrgsolve: mrgsolve models for post-hoc generation of diagnostics and simulation
#model/nonmem/pk: Have pre-run popPK models used in tutorial: Stubs and chains with diagnostic htmls where appropriate

#script
#applied-sims.R : Applied simulations script for pediatric example.
#model-management.R: Generation of models and execution via bbr with commenting on what's being done
#model-diagnostics-report.R: Generation of model diagnostics
#demo-model-table.R: Generation of model tables
#All associated sourced scripts and diagnostics templates etc...


# Read and plot individual parameter estimates from sensitivity an --------

## Read tab files (individual parameter values) 

### Read input data
dat_in <- read_csv(params$data, na=".") 
id <- dat_in %>% select(NUM,ID) %>% distinct()

### Read final model tab file
tab_orig <- map_dfr(seq_len(params$n_chains),function(.chain){
  read_table(file.path(params$modelDir, 
                       params$run, 
                       glue::glue("{params$run}_{.chain}/{params$run}.tab")), 
             skip = 1)}) %>%
  group_by(NUM) %>%
  summarise(across(.fns = median)) %>%
  mutate(sens_num=0) %>% 
  left_join(id)

### Read sensitivity analysis model tab files
tab0 <-  map_dfr(seq_len(params$sens_run), function(.sens_run) {
  
  map_dfr(seq_len(params$n_chains),function(.chain){
    read_table(file.path(params$modelDir, 
                         glue(.sens_run + as.numeric(params$run)), 
                         glue(.sens_run + as.numeric(params$run),"_{.chain}"),
                         glue(.sens_run + as.numeric(params$run),".tab")), 
               skip=1)}) %>%
    group_by(NUM) %>%
    summarise(across(.fns = median)) %>%
    mutate(sens_num=.sens_run)
}) %>% 
  left_join(id)

## remove burn-in iterations and add labels
tab <- tab0 %>%
  bind_rows(tab_orig) %>%
  # select(-NUM) %>%
  distinct() %>%
  mutate(param = ifelse(sens_num == 0, "Final model",      # run sens_num==0 is the final model 
                        "V2"), # runs 1-8 adjusted something affecting V2, 9-16 affected CL
         paramadj = case_when(sens_num == 0 ~ "Final model",
                              sens_num %in% c(2) ~ "Variance 10x",
                              sens_num %in% c(3) ~ "Variance 100x",
                              sens_num %in% c(4) ~ "Increase location 50%",
                              sens_num %in% c(5) ~ "Decrease location 50%", 
                              sens_num %in% c(1) ~ "Fixed"),
         paramadj = factor(paramadj, levels= c("Final model", 
                                               "Fixed",
                                               "Variance 10x", 
                                               "Variance 100x", 
                                               "Increase location 50%",
                                               "Decrease location 50%"
         ))) %>%
  group_by(sens_num) %>%
  mutate(NUM2=1:n()) %>%
  ungroup()

tab2 <- tab %>%
  filter(sens_num==0) %>%
  select(V2orig = V2,ID, NUM2) %>%
  left_join(tab %>% filter(sens_num != 0))

## subset tab data for plotting

### 1). KA variance on V2
tab_ka_var <- tab %>% 
  filter(param %in% c("V2", "Final model")) %>%
  filter(sens_num == 0 | sens_num <= 3)

### 2). KA location on V2
tab_ka_location <- tab %>% 
  filter(param %in% c("V2", "Final model")) %>%
  filter(sens_num == 0 | (sens_num > 3))

## subset tab2 data for plotting

### 1). KA variance on V2
tab2_ka_var <- tab2 %>% 
  filter(param %in% c("V2")) %>%
  filter(sens_num %in% c(2,3))

### 2). KA location on V2
tab2_ka_location <- tab2 %>% 
  filter(param %in% c("V2")) %>%
  filter(sens_num %in% c(4,5))

## Plot 

### 1). plot KA variance on V2 inidividual estimates
tab_ka_var %>%
  ggplot(aes(V2, colour = factor(paramadj), fill = factor(paramadj))) +
  geom_density(alpha=0.1) +
  theme_bw()+
  labs(x = "Individual V2/F Estimate (L)", y = "Density" , colour = "KA prior", fill = "KA prior")

mrggsave_last(stem = "sensitivity-analysis-{params$run}-indV2-KAvar", width = 5, height = 5)

tab2_ka_var %>%
  ggplot() +
  geom_point(aes(x=V2orig,y=V2),alpha=0.5)+
  geom_abline(color="grey",size=1.5)+
  geom_smooth(aes(x=V2orig,y=V2),se=F, formula= y~x,color="blue",linetype=2)+
  facet_wrap(~paramadj)+
  theme_bw()+
  labs(x = "Final Model Individual V2/F Estimate (L)", y = "Sensitivity Analysis Individual V2/F Estimate (L)")

mrggsave_last(stem = "sensitivity-analysis-{params$run}-indV2-KAvar-corr", width = 5, height = 5)

### 2). plot KA location on V2 inidividual estimates
tab_ka_location %>%
  ggplot(aes(V2, colour = factor(paramadj), fill = factor(paramadj))) +
  geom_density(alpha=0.1) +
  theme_bw()+
  labs(x = "Individual V2/F Estimate (L)", y = "Density" , colour = "KA prior", fill = "KA prior")

mrggsave_last(stem = "sensitivity-analysis-{params$run}-indV2-KAscale", width = 5, height = 5)

tab2_ka_location %>%
  ggplot() +
  geom_point(aes(x=V2orig,y=V2),alpha=0.5)+
  geom_abline(color="grey",size=1.5)+
  geom_smooth(aes(x=V2orig,y=V2),se=F, formula= y~x,color="blue",linetype=2)+
  facet_wrap(~paramadj)+
  theme_bw()+
  labs(x = "Final Model Individual V2/F Estimate (L)", y = "Sensitivity Analysis Individual V2/F Estimate (L)")

mrggsave_last(stem = "sensitivity-analysis-{params$run}-indV2-KAscale-corr", width = 5, height = 5)


#Simulate out new PK for atorvastatin example

## Required packages
library(tidyverse)
library(yspec)
library(mrggsave)
library(glue)
library(bbr)
library(here)
library(scales)
library(loo)
library(magrittr)
library(pmtables)
library(patchwork)
library(mrgsolve)

dat_in <- read_csv(here::here("data/derived/atorvWrkShop2.csv"), na=".") 

#Get the duplicate timepoints
olddatdup <- dat_in %>% group_by(ID) %>% mutate(dup = duplicated(TIME)) %>% 
  group_by(ID,TIME) %>% mutate(dup2 = any(dup==T)) %>% ungroup() %>%
  #drop duplicate times
  filter(dup==FALSE)

#load mrgsolve model and simulate design (no LOQ)
modelName <- "2000"
mrgsolve_path = here(glue("model/mrgsolve/{modelName}.mod"))
mod <- mread(mrgsolve_path) %>%
  #Update to 
  param(THETA1 = log(807),
        THETA2 = log(1610),
        THETA3 = log(361),
        THETA4 = log(2030),
        THETA5 = log(0.238)) %>%
  omat(as_bmat(c(0.269,.170,0.895,rep(0,25)))) %>%
  smat(as_dmat(0.150))

#Simulate new data

newdat <- mod %>%
  data_set(olddatdup) %>%
  carry_out(NUM) %>%
  mrgsim(tad=T)

updatedat <- olddatdup %>%
  left_join(newdat %>% select(NUM,Y)) %>%
  mutate(DV = case_when(EVID==0 ~ Y,
                        TRUE ~ DV),
         BLQ = case_when(EVID == 0 & Y < 0 ~ 1,
                         TRUE ~ 0)) %>%
  select(C:NUM)

#Save out the data
data.table::fwrite(
  x = updatedat,
  file = here::here("data", "derived", "atorvWrkShop3.csv"),
  sep = ",",
  quote = FALSE,
  row.names = FALSE,
  na = "."
)








