library(MCMCpack)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(forcats)
library(patchwork)
library(mrggsave)

options(mrggsave.dir = "../deliv/figure", mrg.script = "dfeval.R")
theme_set(pmplots::pm_theme())

scaleMatrix <- c(0.1,0.75)*diag(2)

generate_iw_samples <- function(nsim, .df, .mode) {
  
  sds  <- vector(mode = 'double',length = nsim)
  cors  <- vector(mode = 'double',length = nsim)
  rels  <- vector(mode = 'double',length = nsim)
  vals  <- vector(mode = 'double',length = nsim)
  vals2  <- vector(mode = 'double',length = nsim)
  for (i in 1:nsim) {
    mat <- riwish(.df, .mode)
    sds[i] <- sqrt(mat[1,1])
    cors[i] <- mat[1,2]/sqrt(mat[1,1]*mat[2,2])
    rels[i] <- sqrt(mat[2,2])/.mode[2,2]
    vals[i] <- (mat[1,1])
    vals2[i] <- (mat[1,1])*.df
  }
  data.frame(sd=sds, cor=cors, rel = rels,
             val = vals, val2=vals2)
}

generate_iw_samples2 <- function(nsim, .df, .mode) {
  
  sds  <- vector(mode = 'double',length = nsim)
  cors  <- vector(mode = 'double',length = nsim)
  rels  <- vector(mode = 'double',length = nsim)
  vals  <- vector(mode = 'double',length = nsim)
  vals2  <- vector(mode = 'double',length = nsim)
  for (i in 1:nsim) {
    scaleMatrix1 <- (.df-length(diag(.mode))-1)*.mode
    mat <- riwish(.df, scaleMatrix1)
    sds[i] <- sqrt(mat[1,1])
    cors[i] <- mat[1,2]/sqrt(mat[1,1]*mat[2,2])
    rels[i] <- sqrt(mat[2,2])/.mode[2,2]
    vals[i] <- (mat[1,1])
    vals2[i] <- (mat[2,2])
  }
  data.frame(sd=sds, cor=cors, rel = rels,
             val = vals, val2=vals2)
}

results <- tibble(df=c(4,10,30,100,500))
set.seed(12354)
results1 <- mutate(results, 
                  res = map(df, ~generate_iw_samples(10000, .x, scaleMatrix)))

#Mean with unscaled matrix equal to S/df

results1 %>% 
     mutate(df_c = fct_reorder(paste("df =",df),df)) %>% 
     unnest(res) %>% group_by(df_c) %>% summarise(meanval = mean(val),
                                                  meanval2 = mean(val2))
set.seed(12354)
results2 <- mutate(results, 
                  res = map(df, ~generate_iw_samples2(10000, .x, scaleMatrix)))

results2 %>% 
  mutate(df_c = fct_reorder(paste("df =",df),df)) %>% 
  unnest(res) %>% group_by(df_c) %>% summarise(meanval = mean(val),
                                               meanval2 = mean(val2))
results2 %>% 
  mutate(df_c = fct_reorder(paste("df =",df),df)) %>% 
  unnest(res) %>% 
  pivot_longer(cols=c(sd,cor)) %>% 
  filter(name=='sd') %>% 
  ggplot(aes(x=value, group=df, col=factor(df))) +
  geom_density() +
  facet_wrap(df_c~., scales = 'free') +
  labs(x='Value', y='Density') +
  scale_color_brewer(palette = 'Set1',name='Degrees of freedom') +
  theme_bw()

results2 %>% 
  mutate(df_c = fct_reorder(paste("df =",df),df)) %>% 
  unnest(res) %>% 
  pivot_longer(cols=c(sd,cor,rel)) %>% 
  filter(name=='rel') %>% 
  ggplot(aes(x=value, group=df, col=factor(df))) +
  geom_density() +
  facet_wrap(df_c~., scales = 'free') +
  labs(x='Value', y='Density') +
  scale_color_brewer(palette = 'Set1',name='Degrees of freedom') +
  theme_bw()

p1 <- results2 %>% 
  mutate(df_c = fct_reorder(paste("df =",df),df)) %>% 
  unnest(res) %>% 
  pivot_longer(cols=c(sd,cor,rel,val,val2)) %>% 
  filter(name=='val') %>% 
  ggplot(aes(x=value, group=df, col=factor(df))) +
  geom_density() +
  facet_wrap(df_c~., scales = 'free') +
  labs(x='OMEGA (1,1)', y='Density') +
  scale_color_brewer(palette = 'Set1',name='Degrees of freedom') +
  theme_bw()

p2 <- results2 %>% 
  mutate(df_c = fct_reorder(paste("df =",df),df)) %>% 
  unnest(res) %>% 
  pivot_longer(cols=c(sd,cor,rel,val,val2)) %>% 
  filter(name=='val2') %>% 
  ggplot(aes(x=value, group=df, col=factor(df))) +
  geom_density() +
  facet_wrap(df_c~., scales = 'free') +
  labs(x='Value', y='Density') +
  scale_color_brewer(palette = 'Set1',name='Degrees of freedom') +
  theme_bw()

p3 <- results2 %>% 
  mutate(df_c = fct_reorder(paste("df =",df),df)) %>% 
  unnest(res) %>% 
  pivot_longer(cols=c(sd,cor)) %>% 
  filter(name=='cor') %>% 
  ggplot(aes(x=value, group=df, col=factor(df))) +
  geom_density() +
  labs(x='Correlation (Omega11 ~ Omega22)', y='Density',
       caption = 'A diagonal matrix with diagonal values set at 0.1 was used as the scale matrix') +
  scale_color_brewer(palette = 'Set1',name='Degrees of freedom') +
  theme_bw()


psave1 <- p1 / p3 & theme(legend.position = "bottom") 
psave2 <- psave1 + plot_layout(guides = "collect"); psave2

mrggsave(psave2, tag = "dfinvwisth", height = 8, width=8)
