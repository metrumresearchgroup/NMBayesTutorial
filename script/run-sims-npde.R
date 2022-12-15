#TODO option for point estimate vs full posterior for NPDE/EWRES
#  modelName <- 1100
#  mod_bbr = read_model(file.path(modelDir, modelName, glue("{modelName}_1")))
#  mrgsolve_path = here(glue("script/model/{modelName}.mod"))
#  out_path = file.path(modelDir, modelName, glue("diag-sims-{modelName}.csv"))
#  n_chains = 4; n_post = 1000; ci_level = 95; resid_var = TRUE;
#  point_est = c("median");
#  mat_sqrt = c("cholesky")

#' Generate values for simulation-based diagnostics
#'
#' @param mod_bbr Model to use, either `bbi_nonmem_model` or
#'   `bbi_nonmem_summary` (currently needs to be single chain model, assumed to
#'   be a directory underneath a parent dir, e.g.
#'   "../model/pk/100/100_1/100_1.ctl")
#' @param mrgsolve_path path to mrgsolve model used for simulations
#' @param out_path Path of RDS file to write, relative to current working
#'   directory.
#' @param n_chains Number of chains (default 4)
#' @param n_post Number of posterior samples to use in simulations (default
#'   1000)
#' @param ci_level Level of confidence/credible interval to use (%, default 95)
#' @param resid_var Include residual variability in simulations (default TRUE)
#' @param point_est Type of point estimate to use ("median" (default) or "mean")
#' @param ... Passed to `npde::autonpde()`.
#'   E.g., method for calculating matrix square root: 
#'     * cholesky: decorrelation is performed through the Cholesky decomposition (default)
#'     * inverse: decorrelation is performed by inverting Vi through the eigen function
#'     * polar: the singular-value decomposition (svd) is used
run_sims_npde <- function(
  mod_bbr, mrgsolve_path, out_path,
  n_chains = 4, n_post = 1000, ci_level = 95, resid_var = TRUE,
  point_est = c("median", "mean"), log_dv = FALSE,
  ...
) {
   #browser()
  point_est <- match.arg(point_est)
  
  if (point_est == "median") {
    point_fn <- function(x, ...) median(x, ...)
  } else {
    point_fn <- function(x, ...) mean(x, ...)
  }
  
  qlo <- (1 - ci_level/100) / 2
  qhi <- 1 - qlo
  
  run <- str_remove(get_model_id(mod_bbr), "_\\d+$")
  this_model_dir <- dirname(get_model_path(mod_bbr))
  MODEL_DIR <- dirname(this_model_dir)
  
  tab <- get_chain_files(MODEL_DIR, run, n_chains, "tab",
                         .chain_in_name = FALSE) %>% 
    group_by(NUM) %>% 
    summarise(across(.fns = point_fn)) %>% 
    select(-chain)
  
  data0 <- data.table::fread(get_data_path(mod_bbr), na.strings = ".") %>% 
    as_tibble()
  
  if(isTRUE(log_dv)){
    data0 = data0 %>% mutate(ODV = DV, DV = LDV)
  }
  
  keep <- c("NUM", setdiff(names(tab), names(data0)))
  data <- left_join(tab[,keep], data0, by = "NUM")
  
  
  ## Check no NAs in RATE column
  if("RATE" %in% names(data)){ 
    data <- data %>% mutate(RATE = ifelse(is.na(RATE), 0, RATE))
  }
  ## Check no negatives in SS column
  if("SS" %in% names(data)){ 
    data <- data %>% mutate(SS = ifelse(-999, 0, SS))
  }
  ## Check if an EVID column exists, if not, create one
  if(!("EVID" %in% names(data))){ 
    data <- data %>% mutate(EVID = 0)
  }
  
  
  ext <- get_chain_files(MODEL_DIR, run, n_chains, "ext") %>% 
    filter(ITERATION > 0)
  
  require("mrgsolve")
  mod_mrgsolve <- mread(mrgsolve_path)
  mod_mrgsolve <- update(mod_mrgsolve, rtol = 1e-5, atol = 1e-5, outvars = "Y")
  
  if (!resid_var) mod_mrgsolve <- mod_mrgsolve %>% mrgsolve::zero_re("sigma")
  
  # simulations using
  #   - random samples from posterior estimates
  #   - random samples from between- and within-subject variability distributions
  set.seed(1)
  post <- ext %>% 
    slice_sample(n = n_post) %>% 
    select(-MCMCOBJ) %>% 
    mutate(sample = row_number())
  
  # EPRED -------------------------------------------------------------------
  
  sim_pred <- furrr::future_map_dfr(
    post$sample,
    .options = opt,
    .progress = interactive(),
    .f = function(.sample) {
      output <- mod_mrgsolve %>% 
        mrgsolve::carry_out(NUM) %>% 
        ##  intention here is not to include the post hoc ETAs, 
        ##  only ETAs simulated from OMEGA (below)
        mrgsolve::data_set(data %>% select(-contains("ETA"))) %>%
        mrgsolve::param(
          post %>% 
            filter(sample == .sample) %>%
            select(starts_with("THETA"))
        ) %>%
        mrgsolve::omat(mrgsolve::as_bmat(post %>% filter(sample == .sample), "OMEGA")) %>%
        ### Needs to be as_bmat not as_dmat
        mrgsolve::smat(mrgsolve::as_bmat(post %>% filter(sample == .sample), "SIGMA")) %>%
        mrgsolve::mrgsim_df(obsonly = TRUE) 
      
      output <- output %>% 
        select(NUM, DV_sim = Y) %>% 
        mutate(sample = .sample)
 
      return(output)
    }
  ) %>% 
    as_tibble()
  
  data <- data %>% 
    select(-any_of("EPRED")) %>% 
    left_join(
      sim_pred %>% 
        group_by(NUM) %>% 
        summarise(
          EPRED = point_fn(DV_sim),
          EPRED_lo = quantile(DV_sim, qlo),
          EPRED_hi = quantile(DV_sim, qhi)
        ),
      by = "NUM"
    )
  
  # IPRED -------------------------------------------------------------------
  
  #TODO modify build_path_from_model to work with chains
  ipred_sim <- FALSE
  fn_ipar <- file.path(
    MODEL_DIR,
    run,
    glue::glue("{run}_{n_chains}"),
    glue::glue("{run}_{n_chains}.iph")
  )
  if (file.exists(fn_ipar)) {
    ipred_sim <- TRUE
    ipar <- get_chain_files(MODEL_DIR, run, n_chains, "iph", verbose = TRUE) %>% 
      filter(ITERATION > 0) %>% 
      select(chain, ITERATION, ID, starts_with("ETA"))
    
    # add median ETAs to table
    data <- data %>% 
      select(-starts_with("ETA")) %>% 
      left_join(
        ipar %>% 
          select(-chain, -ITERATION) %>% 
          # rename from "ETA(1)" to "ETA1" etc
          rename_with(~ str_replace_all(.x, "[()]", ""), starts_with("ETA")) %>% 
          group_by(ID) %>% 
          summarise(across(starts_with("ETA"), point_fn)),
        by = "ID"
      )
    
    # simulations using
    #   - random samples from individual posterior estimates (post hoc ETAs)
    set.seed(1)
    post_ipar <- ipar %>% 
      # add corresponding THETA and SIGMA estimates for each posterior sample
      left_join(
        post %>% select(-starts_with("OMEGA")),
        by = c("chain", "ITERATION")
      ) %>% 
      filter(!is.na(sample)) %>% 
      select(-chain) %>% 
      # rename from "ETA(1)" to "ETA1" etc
      rename_with(~ str_replace_all(.x, "[()]", ""), starts_with("ETA"))
    print("Simulate IPREDs")
    
    sim_ipred <- furrr::future_map_dfr(
      unique(post_ipar$sample),
      .options = opt,
      .progress = interactive(),
      .f = function(.sample) {
        output <- mod_mrgsolve %>% 
          mrgsolve::zero_re("omega") %>% 
          mrgsolve::carry_out(NUM) %>% 
          mrgsolve::data_set(data) %>% 
          mrgsolve::idata_set(post_ipar %>% filter(sample == .sample)) %>% 
          ### Needs to be as_bmat not as_dmat
          mrgsolve::smat(mrgsolve::as_bmat(post %>% filter(sample == .sample), "SIGMA")) %>%
          mrgsolve::mrgsim_df(obsonly = TRUE) 
        
        output <- output %>% 
          select(NUM, DV_sim = Y) %>% 
          mutate(sample = .sample)
        
#       Moved back-transformation code below - caused bias in residuals to back-transform here
        # if(isTRUE(log_dv)){ 
        #   output <- output %>%
        #     mutate(DV_sim = exp(Y)) %>% 
        #     select(NUM, DV_sim) %>% 
        #     mutate(sample = .sample)
        # } else {
        #   output <- output %>% 
        #     select(NUM, DV_sim = Y) %>% 
        #     mutate(sample = .sample)
        # }
        
        return(output)
      }
    ) %>% 
      as_tibble()
    
    data <- data %>% 
      select(-any_of("IPRED")) %>% 
      left_join(
        sim_ipred %>% 
          group_by(NUM) %>% 
          summarise(
            IPRED = point_fn(DV_sim),
            IPRED_lo = quantile(DV_sim, qlo),
            IPRED_hi = quantile(DV_sim, qhi)
          ),
        by = "NUM"
      )
    
  }
  
  # EWRES and NPDE ----------------------------------------------------------
  
  df_obs <- data %>%
    filter(EVID == 0) %>%
    select(ID, TIME, DV, NUM)
  
  df_sim <- sim_pred %>%
    left_join(
      data %>% filter(EVID == 0) %>% select(NUM, ID, TIME),
      by = "NUM"
    ) %>% 
    select(ID, TIME, DV = DV_sim)
  
  fn_df_obs <- file.path(tempdir(), "df_obs.txt")
  fn_df_sim <- file.path(tempdir(), "df_sim.txt")
  write_delim(df_obs, fn_df_obs)
  write_delim(df_sim, fn_df_sim)
  npde_out <- npde::autonpde(
    namobs = fn_df_obs,
    namsim = fn_df_sim,
    iid = "ID", ix = "TIME", iy = "DV",
    boolsave = FALSE,
    ...
  )
  ewres_npde <- npde_out@results@res %>% 
    bind_cols(df_obs %>% select(NUM)) %>% 
    select(NUM, EWRES = ydobs, NPDE = npde)
  
  data <- data %>% 
    select(-any_of(c("EWRES", "NPDE"))) %>% 
    left_join(ewres_npde, by = "NUM")
  
  keep <- c("NUM", setdiff(names(data), names(data0)))
  data <- data[, keep]
  

  # Back transform if model estimated in the log domain
  if(isTRUE(log_dv)){
    data <- mutate(data, 
                    
                    LNIPRED = IPRED, 
                    LNIPRED_lo = IPRED_lo, 
                    LNIPRED_hi = IPRED_hi, 
                    LNEPRED = EPRED,
                    LNEPRED_lo = EPRED_lo,
                    LNEPRED_hi = EPRED_hi,
                    
                    IPRED = exp(IPRED),
                    IPRED_lo = exp(IPRED_lo),
                    IPRED_hi = exp(IPRED_hi),
                    EPRED = exp(EPRED),
                    EPRED_lo = exp(EPRED_lo),
                    EPRED_hi = exp(EPRED_hi)
    )  
  }
  
  saveRDS(
    list(
      data = data,
      mod_bbr = mod_bbr,
      mrgsolve_path = mrgsolve_path,
      out_path = out_path,
      n_chains = n_chains,
      n_post = n_post,
      ci_level = ci_level,
      resid_var = resid_var,
      point_est = point_est,
      ipred_sim = ipred_sim
    ),
    file = out_path
  )
}
