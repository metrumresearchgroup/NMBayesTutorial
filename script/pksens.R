# Packages, params and options --------------------------------------------

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

runno <- '2000'

params <- list(
  run = runno, 
  n_chains = 4,
  sens_run = 5, # number of sensitivity analysis runs
  script = "pksens.R",
  modelDir = here::here("model/pk"),
  tables = here::here("deliv/table/"), 
  figures = here::here("deliv/figure/"),
  #spec = "../script/analysis.yml", 
  data = here::here("data/derived/atorvWrkShop3.csv")
)

if (!dir.exists(params$tables)) dir.create(params$tables)
if (!dir.exists(params$figures)) dir.create(params$figures)

figdir_run <- file.path(params$figures, runno)
if (!dir.exists(figdir_run)) dir.create(figdir_run) 
tabdir_run <- file.path(params$tables, 'pk', runno)
if (!dir.exists(tabdir_run)) dir.create(tabdir_run) 

options(mrggsave.dir =  figdir_run)
options(mrg.script = params$script)


# Helper functions --------------------------------------------------------
source(here::here("script/functions-diagnostics.R"))


# Read and plot posterior estimates from sensitivity analyses -------------

## Read ext files of the final model (parameter values)
ext_orig <- map_dfr(seq_len(params$n_chain), function(.chain) {
  data.table::fread(
    file = file.path(params$modelDir, params$run, glue("{params$run}_{.chain}"),
                     glue("{params$run}_{.chain}.ext"))
  ) %>% 
    as_tibble() %>% 
    mutate(chain = .chain)
}) %>%
  mutate(sens_num = 0)

## Read ext files of the sensitivity analysis model
ext0 <- map_dfr(seq_len(params$sens_run), function(.sens_run) {
  
  map_dfr(seq_len(params$n_chain), function(.chain) {
    data.table::fread(
      file = file.path(params$modelDir, glue(.sens_run + as.numeric(params$run)), 
                       glue(.sens_run + as.numeric(params$run),"_{.chain}"),
                       glue(.sens_run + as.numeric(params$run),"_{.chain}.ext"))
    ) %>% 
      as_tibble() %>% 
      mutate(chain = .chain)
  })  %>% 
    
    mutate(sens_num = .sens_run)
})

#One with KA fixed (2001)
#One with KA prior inflated 10x (2002)
#One with KA prior variance 100x (2003)
#One with KA prior mean increased 50% (2004)
#One with KA prior mean decreased 50% (2005)
## remove burn-in iterations and add labels
ext <- ext0 %>%  
  bind_rows(ext_orig) %>%
  filter(ITERATION > 0) %>%
  mutate(param = ifelse(sens_num == 0, "Final model",     # run sens_num==0 is the final model 
                        "V2"), # runs 1-5 adjusted something affecting V2
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
                                               )))


## subset data for plotting

### 1). KA variance on V2
ext_ka_var <- ext %>% 
  filter(param %in% c("V2", "Final model")) %>%
  filter(sens_num == 0 | sens_num <= 3) %>%
  gather(THETA1:THETA2,key='key',value = 'value') %>%
  mutate(key = case_when(key=='THETA1' ~ 'CL/F (L/h)',
                         key=='THETA2' ~ 'V2/F (L)',
                         TRUE ~ key))

### 2). KA location on V2
ext_ka_location <- ext %>% 
  filter(param %in% c("V2", "Final model")) %>%
  filter(sens_num == 0 | (sens_num > 3)) %>%
  gather(THETA1:THETA2,key='key',value = 'value') %>%
  mutate(key = case_when(key=='THETA1' ~ 'CL/F (L/h)',
                         key=='THETA2' ~ 'V2/F (L)',
                         TRUE ~ key))

## Plot 

### 1). plot KA variance on density
p1 <- ext_ka_var %>%
  ggplot(aes(exp(value), colour = factor(paramadj), fill = factor(paramadj))) +
  geom_density(alpha=0.1) +
  theme_bw()+
  labs(x = "Typical Estimate", y = "Density" , 
       colour = "KA prior", fill = "KA prior") +
  facet_wrap(~key,scales = 'free')

### 2). plot KA scale on V2 density
p2 <- ext_ka_location %>% 
  ggplot(aes(exp(value), colour = factor(paramadj), fill = factor(paramadj))) +
  geom_density(alpha=0.1) +
  theme_bw()+
  labs(x = "Typical Estimate", y = "Density" , 
       colour = "KA prior", fill = "KA prior") +
  facet_wrap(~key,scales = 'free')


psave <- p1/p2; psave

mrggsave_last(stem = "sensitivity-analysis-{params$run}-V2-CL-KA", width = 7, height = 7)



# Calculate PSS-LOO for each sensitivity analysis run ---------------------

## original chain  

iph <- get_chain_files(params$modelDir, params$run, params$n_chains, "iph") %>% 
  filter(ITERATION > 0)

iobj <- iph %>% select(chain,ITERATION, ID, MCMCOBJ)

n_iter <- max(iobj$ITERATION)
n_id <- n_distinct(iobj$ID)
iobj_array <- array(
  double(n_iter * params$n_chains * n_id),
  dim = c(n_iter, params$n_chains, n_id)
)
for (.chain in seq_len(params$n_chains)) {
  iobj_array[,.chain,] <- iobj %>% 
    filter(chain == .chain) %>% 
    arrange(ID, ITERATION) %>% 
    pull(MCMCOBJ) %>% 
    as.matrix()
}

rel_n_eff <- relative_eff(exp(iobj_array))
loo_orig <- loo(iobj_array, r_eff = rel_n_eff) 

# loo_orig<- loo_out$estimates %>% as_tibble() %>% slice(3) %>% mutate(sens_num=0)

#store ELPD of final model
loo_orig_elpd <- loo_orig$estimates %>% as_tibble(rownames="metric")  %>% filter(metric=="elpd_loo") 
loo_orig_elpd <- loo_orig_elpd$Estimate

## sensitivity analysis runs
loo_sens <- 
  map_dfr(seq_len(params$sens_run), function(.sens_run) {
    

    iph <- get_chain_files(params$modelDir,
                           glue(.sens_run + as.numeric(params$run)),
                           params$n_chains, "iph") %>% 
      filter(ITERATION > 0)
    
    iobj <- iph %>% select(chain, ITERATION, ID, MCMCOBJ)
    
    n_iter <- max(iobj$ITERATION)
    n_id <- n_distinct(iobj$ID)
    iobj_array <- array(
      double(n_iter * params$n_chains * n_id),
      dim = c(n_iter, params$n_chains, n_id)
    )
    for (.chain in seq_len(params$n_chains)) {
      iobj_array[,.chain,] <- iobj %>% 
        filter(chain == .chain) %>% 
        arrange(ID, ITERATION) %>% 
        pull(MCMCOBJ) %>% 
        as.matrix()
    }
    
    rel_n_eff <- relative_eff(exp(iobj_array))
    loo_out <- loo(iobj_array, r_eff = rel_n_eff) #looic is the LOO information criterion
    
    # loo_out$estimates %>% as_tibble() %>% slice(3) %>% mutate(sens_num=.sens_run)
    print(loo_out)
    
    x <- print(loo_compare(loo_out, loo_orig), simplify=F, digits=3) %>% 
      as_tibble() %>% 
      mutate(across(everything(),~as.numeric(.x)),
             sens_num=.sens_run)
    
    x
  })

# formatting so final model is the reference model
# for each pairwise comparison, 
# if the max value equals the final model, retain the sign of elpd_diff, otherwise multiply it by -1
orig_row <- loo_sens %>% 
  filter(sens_num == 1 & elpd_loo == loo_orig_elpd) %>% 
  mutate(sens_num = 0) %>%
  mutate(elpd_diff = 0,
         se_diff = 0)

comp<- loo_sens %>%
  mutate(elpd_diff = ifelse(elpd_loo!=loo_orig_elpd,elpd_diff,elpd_diff*-1),
         elpd_loo = ifelse(elpd_loo==loo_orig_elpd,NA_real_,elpd_loo),
         se_elpd_loo = ifelse(elpd_loo==loo_orig_elpd,NA_real_,se_elpd_loo),) %>%
  group_by(sens_num) %>%
  fill(c(elpd_loo,se_elpd_loo),.direction = c("updown")) %>%
  ungroup() %>%
  filter(elpd_diff != 0) %>%
  bind_rows(orig_row) %>%
  arrange(sens_num) %>%
  mutate(across(c(-sens_num), ~sig(.x,digits = 4))) %>%
  select(sens_num, elpd_loo, se_elpd_loo, elpd_diff, se_diff)

loo_all <- comp %>%
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
         )))

## Output table for loo evaluation

looTab = loo_all %>% 
  mutate(group = case_when(sens_num == 0 ~ "",
                           param == "V2"&sens_num %in% c(1,2,3,4,5) ~ "KA prior adjustment")) %>%
  mutate(group = fct_relevel(group, c("", 
                                      "KA prior adjustment"))) %>% 
  arrange(group) %>% 
  select(-sens_num,-param) %>%
  relocate(paramadj) %>% 
  stable(
    cols_rename=vars("ELPD"="elpd_loo",
                     "Standard error" ="se_elpd_loo",
                     "Difference" = "elpd_diff",
                     "Standard error ... of difference" = "se_diff",
                     "Model" = "paramadj"),
    hline_at=2,
    panel = "group", 
    clear_reps="group",
    r_file=params$script,
    output_file=file.path(params$tables, "LOO-summary.tex"),
    notes = "Abbreviations: ELPD: expected log predictive density",
    align = cols_center(group = col_ragged(8)),
    lt_cap_text = 'Sensitivity Analysis: Summary of ELPD values',
    lt_cap_label = "tab:pk-loo-cv"
  ) 

# looTab %T>% st2report(ntex=2) 
looTab %T>% stable_save(dir = tabdir_run, file = "LOO-summary.tex")

st2doc(
  list(looTab),
  output_dir = tabdir_run,
  output_file = glue("LOO-summary-preview.pdf")
)

