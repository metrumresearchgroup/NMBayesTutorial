### General table-related helper functions ####  -------------------------------------------------------

asNum  <- function(f){ return(as.numeric(as.character(f))) }
sig = pmtables::sig
select <- dplyr::select
rename <- dplyr::rename
ggplot2::theme_set(ggplot2::theme_bw())
parens <- function(x) paste0('(',x,')')
parensSQ <- function(x) paste0('[',x,']')
parensSQ_CV <- function(.x) glue::glue("[CV\\%=<<.x>>]", .open = "<<", .close  = ">>")
parensSQ_corr <- function(.x) glue::glue("[Corr=<<.x>>]", .open = "<<", .close  = ">>")
parensSQ_se <- function(.x) glue::glue("[SD=<<.x>>]", .open = "<<", .close  = ">>")
getEvenNo = function(x) x[which(x %% 2 == 0)]


# when using bootstrap estimates
qt   <- function(x, prob) unname(quantile(x, probs = prob, na.rm = TRUE))
med  <- function(x) qt(x, 0.5)
lo  <- function(x) qt(x, 0.025)
hi  <- function(x) qt(x, 0.975)


### For parameter tables ####  -------------------------------------------------------
# CV equations https://ascpt.onlinelibrary.wiley.com/doi/full/10.1002/psp4.12404
getCV_propO  <- function(v) sqrt(exp(v) - 1) * 100
getCV_propS  <- function(v) sqrt(v) * 100
lowerCI <- function(est, se) est - 1.96*se
upperCI <- function(est, se) est + 1.96*se
parensSQ <- function(x) paste0('[',x,']')

mathMode <- function(.x) glue::glue("$<<.x>>$", .open = "<<", .close  = ">>")
gtGreek <- function(.x) glue::glue("\\<<.x>>", .open = "<<", .close  = ">>")
greekNum <- function(.x, .y) glue::glue("<<.x>>_{<<.y>>}", .open = "<<", .close  = ">>")
expGreek  <- function(.x, .y) glue::glue("$\\exp(\\<<.x>>_{<<.y>>})$", .open = "<<", .close  = ">>")



## Define a series of true/false columns to make filter easier later
defineRows <- function(df){
  df %>% 
    mutate(
      TH = ifelse(stringr::str_detect(name, "TH"), T, F),
      OM = ifelse(stringr::str_detect(name, "OM"), T, F),
      S = ifelse(stringr::str_detect(name, "S"), T, F),
      LOG = ifelse(name %in% logParam, T, F),
      propErr = ifelse(name %in% propErr, T, F),
      addErr = ifelse(name %in% addErr, T, F))
}


## calculate 95% confidence intervals
get95CI <- function(df){
  df %>%      
    mutate(lower = lowerCI(value, se), 
           upper = upperCI(value, se))
}

## calculate % RSE - not used but included if needed 
# Note, this is appropriate when parameters are estimated untransformed or in the log 
# it may not be appropriate if any other transformations (such as logit) were used
getpRSE <- function(df){
    df %>%
      mutate(pRSE = case_when(TH & LOG ~ sig ((exp(se^2)-1)*100), # suggested by JF
                              TH  ~ sig ((se/abs(value)) * 100),
                              diag ~ sig ((se/abs(value)) * 100),
                              TRUE ~ "-"))
}

## Back transform parameters estimated in the log domain
# make sure any other calculations, such as CI (and pRSE) are 
# done before back-calculating these values
backTrans_log <- function(df){
  df %>% 
    mutate(value = case_when(LOG ~ exp(value), TRUE ~ value),
           lower = case_when(LOG ~ exp(lower), TRUE ~ lower),
           upper = case_when(LOG ~ exp(upper), TRUE ~ upper))
}

# value should have estimate [something]
#   theta = estimate only                           # use estimate column
#   omega diagonals = variance [%CV]                # estimate [CV from estimate, stderr]
#   omega off-diagonals = covariance [corr coeff]   # estimate [random_effect_sd]
#   sigma diagonal proportional = variance [%CV]    # estimate [CV from estimate, stderr]
#   sigma diagonal additive = variance [SD]         # estimate [random_effect_sd]
getValueSE <- function(df){
  df %>% 
    mutate(value = estimate, 
           se = stderr,
           corr_SD = case_when(OM & !diag |
                                 S & diag & addErr ~ sig(random_effect_sd), 
                               TRUE ~ "-")
    )
}


# 95% CI should show lower, upper or FIXED
# rounding for display in report
# define what is in estimate column and what is in square brackets 
formatValues <- function(df){
  df %>% 
    mutate(ci = paste0(sig(lower), ', ', sig(upper)), 
           ci = if_else(fixed, "FIXED", ci),
           # get % CV 
           cv = case_when(diag & OM ~ sig(getCV_propO(value)), 
                          diag & S & propErr ~ sig(getCV_propS(value)), 
                          TRUE ~ "-"),
           # round values for report table
           value = sig(value),           
           
           # define which values appear where
           value = case_when(diag & OM | 
                               diag & S & propErr ~ 
                               glue::glue("{value} {parensSQ_CV(cv)}"),
                             !diag & OM ~ glue::glue("{value} {parensSQ_corr(corr_SD)}"),
                             diag & S & addErr ~ glue::glue("{value} {parensSQ_se(corr_SD)}"),
                             !diag & S ~ glue::glue("{value} {parensSQ_corr(corr_SD)}"),
                             TRUE ~ value),
           # round shrinkage values for report table
           shrinkage = case_when(is.na(shrinkage) ~ "-", 
                              TRUE ~ sig(shrinkage)))
}

## Format the THETA/OMEGA/SIGMA values to display as greek letters with 
# subscript numbers
formatGreekNames_2.0 <- function(df){
  df %>% 
    mutate(greekName = name) %>% 
    # make column with greek letters and parameter numbers
    separate(greekName,
             into = c("text", "num"),
             sep = "(?<=[A-Za-z])(?=[0-9])"
    ) %>% 
    separate(parameter_names,
             into = c("text2", "num2"),
             sep = "A"
    ) %>% 
    select(-num, -text2) %>% 
    mutate(text = case_when(OM ~ "Omega", 
                            S ~ "Sigma",
                            TRUE ~ tolower(text)),
           greek = case_when(TH & LOG ~ expGreek(text, num2),
                             TRUE ~ mathMode(greekNum(gtGreek(text), num2))
           )
    )
}



## Define which parameters should appear under which panel name in the final table
getPanelName = function(df){
  rbind(df %>% filter(TH), df %>% filter(OM), df %>% filter(S)) %>% 
    mutate(type = case_when(propErr | addErr ~ "Residual variance",
                            OM & !diag ~ "Interindividual covariance parameters",
                            str_detect(abb, "IIV") ~ "Interindividual variance parameters",
                            # IOV not used here but included for convenience 
                            str_detect(abb, "IOV") ~ "Interoccasion variance parameters",  
                            abb %in% covarList ~ "Covariate effect parameters",
                            str_detect(name, "TH") & !(name %in% covarParams) ~ "Structural model parameters"),
           # Make type a factor and use to sort, this ensures all parameters 
           # of the same type are together - needed to make sure panels pull out 
           # correct rows
           type_f = factor(type, levels=unique(type))
    ) %>% 
    arrange(type_f)
}


# Format values for bootstrap run
formatValues_boot <- function(df){
  df %>% 
    mutate(boot_ci = paste0(sig(lower), ', ', sig(upper)), 
           boot_ci = if_else(fixed, "FIXED", boot_ci),
           boot_value = sig(value))
}


## If not using bbr function: extract required parameters from the ext file ---------------------------
read_extfile <- function(file) {
  tab <- read_table(file,skip = 1)  
  iteration <- tab[["ITERATION"]]
  tab[c("ITERATION", "OBJ")] <- NULL
  estim <- filter(tab, iteration == -1000000000)        # Estimate
  stderr <- filter(tab,    iteration == -1000000001)    # Standard error of the estimate
  random_effect_sd <- filter(tab,    iteration == -1000000004)  # Correlation matrix
  random_effect_sdse<- filter(tab,    iteration == -1000000005) # SE of correlation matrix
  tibble(
    parameter_names = names(estim), 
    estimate = unlist(estim),
    stderr = unlist(stderr),
    random_effect_sd = unlist(random_effect_sd),
    random_effect_sdse = unlist(random_effect_sdse),
    diag = str_match_all(parameter_names, "([0-9]+),([0-9]+)") %>%
      map_lgl(~.x[2]==.x[3]) %>% replace_na(FALSE),
    fixed = stderr > 1e9
  )
}

## If not using bbr function: extract shrinkage from the shk file ---------------------------
read_shkfile <- function(file) {
  tab <- suppressWarnings(read_table(file,skip = 1))  
  shrink <- tab %>% 
    filter(X2 == 4) %>% 
    select(names(tab)[str_detect(colnames(tab), "ET")])        
  tibble(
    etas = names(shrink), 
    num = unlist(stringr::str_extract(etas, "\\d+")),
    name = paste0("OMEGA", num, num),
    shrinkage = unlist(shrink),
  ) %>% 
    select(name, shrinkage)
}




## For study summary

beginItemized <- function() return("\\begin{itemize}")
endItemized   <- function() return("\\end{itemize}")
item          <- function() return("\\item")
bold          <- function(.x) glue("\\textbf{<<.x>>}", .open = "<<", .close  = ">>")
contents      <- function(.x) paste0(names(.x), ": ", .x)

perStudy = function(...){
  row_list = list(...)
  
  # grab file extension from first element and then remove it from the list
  fileEx = row_list[1]
  row_list = row_list[-1]
  
  # First element is the upper level of the list
  out = c(beginItemized(),
          item(), bold(contents(row_list[1])),
          beginItemized())
  
  # make second level of bullets
  for(rr in 2:length(row_list)){
    if(!is.na(row_list[rr])) {
      out = c(out, item(), contents(row_list[rr]))
    }
  }
  
  # end both itemized lists
  out = c(out, endItemized(),
          endItemized()
  )
  
  # save the information out as characters
  out %>%
    as.character() %>%
    writeLines(con = file.path(tabDir, glue("textStudySum_{fileEx}.tex")))
}


