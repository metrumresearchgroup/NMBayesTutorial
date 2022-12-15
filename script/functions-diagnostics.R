##' Purpose: 
##' define functions that are called repeatedly while producing diagnostics
##' for model runs
##' 

##' map_wrap_eta_cont
##' Purpose: plot ETA vs all continuous covariates
##' @param  .map_etas character: name of ETA
##' @param  .co       string: continuous covariate list
##' @param  .id       dataframe: includes one record per id
map_wrap_eta_cont <- function(.map_etas,.co,.id) {
  .p <- wrap_eta_cont(
    .id,
    y = .map_etas,
    x = .co,
    use_labels = TRUE,
    ncol = 2, scales = "free_x"
  )
}

##' map_eta_cat
##' Purpose: plot all ETAs vs a categorical covariate
##' @param  .map_ca   character: name of a categorical covariate
##' @param  .etas     string: ETA list
##' @param  .id       dataframe: includes one record per id
map_eta_cat <- function(.map_ca, .etas, .id) {
  .p <- eta_cat(.id, x = .map_ca, y = .etas) %>% 
    ## CHECK: depending on the labels, this may need to be changed 
    purrr::map(~.x+rot_x(45)) %>% 
    pm_grid
}


trace_plot <- function(.ext_plot, param_per_plot = 9, nrow = 3, ncol = 3,
                       alpha = 1) {
  ans <- pmplots:::chunk_by_id(.ext_plot, nchunk = param_per_plot,
                               id_col = "name")
  p <- lapply(ans, function(.x) {
    .x %>% 
      ggplot(aes(ITERATION, value, colour = chain)) +
      geom_line(alpha = alpha) +
      facet_wrap(~ name, scales = "free_y", nrow = nrow, ncol = ncol) +
      labs(
        x = "Iteration",
        y = "Value",
        colour = "Chain"
      )
  })
  return(p)
}

density_plot <- function(.ext_plot, param_per_plot = 9, nrow = 3, ncol = 3,
                         by_chain = TRUE) {
  ans <- pmplots:::chunk_by_id(.ext_plot, nchunk = param_per_plot,
                               id_col = "name")
  p <- lapply(ans, function(.x) {
    if (by_chain) {
      tmp <- .x %>% ggplot(aes(value, colour = chain)) +
        labs(
          x = "Value",
          y = "Density",
          colour = "Chain"
        ) +
        geom_density()
    } else {
      tmp <- .x %>% ggplot(aes(value)) +
        labs(
          x = "Value",
          y = "Density"
        ) +
        geom_density(fill = "lightblue")
    }
    tmp <- tmp +
      geom_density() +
      facet_wrap(~ name, scales = "free", nrow = nrow, ncol = ncol)
    return(tmp)
  })
  return(p)
}

#' Calculate shrinkage (baysh): Calculate shrinkage of posterior means
#' @param .ipar, dataframe, posterior of ETA samples
#' @param .ext, dataframe, posterior of OMEGA (and other) estimates
#' @export
baysh <- function(.ipar, .ext) {
  #Gelman & Pardoe
  num <- .ipar %>% 
    group_by(ID) %>% 
    select(ID, starts_with("ETA")) %>% 
    summarise(across(everything(), mean), .groups = "drop") %>% 
    select(starts_with("ETA")) %>% 
    summarise(across(everything(), sd))
  eta_nums <- str_extract(names(num), "\\d+")
  omegas <- map_chr(eta_nums, ~ glue::glue("OMEGA({.x},{.x})"))
  denom1 <- .ext %>% 
    select(all_of(omegas)) %>% 
    summarise(across(everything(), ~ sqrt(mean(.x))))
  shrinkage <- as.numeric(1 - num/denom1) * 100
  out <- tibble(eta = names(num), shrinkage = shrinkage)
  return(out)
}

#' Get output files from multiple Bayes chains
#'
#' @param .model_dir high-level model directory, e.g. "../../model/pk"
#' @param .run run number
#' @param .n_chains number of chains
#' @param .ext file extension to read (e.g., "ext" or "tab")
#' @param .chain_in_name logical, TRUE (default) if file name includes the chain
#'   number (e.g., "1000_1.ext" vs "1000.tab")
#' @param .chain_delim delimiter used to separate chain number in file name
#'   (default "_")
#'
#' @return Tibble including all rows from each chain output file, with
#'   additional column `chain` indicating chain number
#' @export
#'
#' @examples
#'   get_chain_files(MODEL_DIR, params$run, params$n_chains, "ext")
get_chain_files <- function(.model_dir, .run, .n_chains, .ext,
                            .chain_in_name = TRUE, .chain_delim = "_",
                            verbose = FALSE) {
#TODO check that these files exist first
  map_dfr(seq_len(.n_chains), function(.chain) {
    if (verbose) print(glue::glue("Reading chain {.chain}..."))
    if (.chain_in_name) {
      fn <- glue::glue("{.run}{.chain_delim}{.chain}.{.ext}")
    } else {
      fn <- glue::glue("{.run}.{.ext}")
    }
    data.table::fread(
      file = file.path(
        .model_dir,
        .run,
        glue::glue("{.run}{.chain_delim}{.chain}"),
        fn
      )
    ) %>% 
      as_tibble() %>% 
      mutate(chain = .chain)
  })
}

sqrt_chol <- function(x) {
  sqrt_mat <- t(chol(x))
  return(sqrt_mat)
}

sqrt_eigen <- function(x) {
	var.eig <- eigen(x)
	sqrt_mat <- var.eig$vectors %*%
	  diag(sqrt(var.eig$values)) %*%
	  solve(var.eig$vectors)
  return(sqrt_mat)
}

sqrt_polar <- function(x) {
	xmat <- chol(x)
	svdec <- svd(xmat)
	umat <- svdec$u %*% t(svdec$v)
	sqrt_mat <- t(umat) %*% xmat
	return(sqrt_mat)
}

#' Calculate EWRES and NPDE
#'
#' @param .df_obs observed data frame containing:
#' * `NUM` (row number in original dataset)
#' * `ID`
#' * `DV`
#' @param .df_sim simulated data frame containing:
#' * `NUM`
#' * `DV_sim`
#' * `sample` (sample or replicate number)
#' @param .sqrt method for calculating matrix square root:
#' * eigen: eigenvalue decomposition (default)
#' * chol: Cholesky decomposition
#' * polar: polar decomposition
#'
#' @return data frame containing columns `EWRES`, `NPDE`, and `NUM`.
#' @export
#'
#' @examples
calc_ewres_npde <- function(.df_obs, .df_sim,
                            .sqrt = c("eigen", "chol", "polar")) {
  #TODO npde for BQL
  .sqrt <- match.arg(.sqrt)
  .df <- left_join(.df_obs, .df_sim, by = "NUM")
  n_sample <- n_distinct(.df_sim$sample)
  
  .out <- furrr::future_map_dfr(
    unique(.df$ID),
    .options = opt,
    .progress = interactive(),
    function(.id) {
      tmp <- filter(.df, ID == .id)
      tmp0 <- filter(tmp, sample == first(tmp$sample))
      # column for each observation within subject
      # row for each sample
      matsim <- matrix(tmp$DV_sim, nrow = n_sample)
      evar <- cov(matsim)
      #TODO change to median?
      #    epred <- matsim %>%
      #      as_tibble() %>%
      #      summarise(across(.fns = median)) %>% 
      #      pivot_longer(cols = everything()) %>% 
      #      pull(value)
      epred <- colMeans(matsim)
      if (det(evar) <= 0 | any(eigen(evar)$values <= 0)) {
        y_star <- NA
        npde <- NA
      } else {
        if (.sqrt == "chol")  evar_sqrt <- sqrt_chol(evar)
        if (.sqrt == "eigen") evar_sqrt <- sqrt_eigen(evar)
        if (.sqrt == "polar") evar_sqrt <- sqrt_polar(evar)
        evar_sqrt_inv <- solve(evar_sqrt)
        y_star <- evar_sqrt_inv %*% (tmp0$DV - epred)
        ysim_star <- map(unique(tmp$sample), function(.sample) {
          #tmp_k <- filter(tmp, sample == .sample)
          # this is > 10x faster than filter! why?
          tmp_k <- tmp[tmp$sample == .sample, ]
          ysim_star_k <- evar_sqrt_inv %*% (tmp_k$DV_sim - epred)
          return(ysim_star_k)
        }) %>% 
          unlist()
        
        if (is.complex(y_star)) {
          y_star <- NA
        }
        
        if (is.complex(ysim_star)) {
          npde <- NA
        } else {
          delta_star <- matrix(ysim_star < rep(y_star, n_sample), ncol = n_sample)
          npde <- qnorm(rowMeans(delta_star))
        }
      }
      .out_id <- tibble(
        EWRES = as.numeric(y_star), 
        NPDE  = as.numeric(npde),
        NUM   = tmp0$NUM
      )
      return(.out_id)
    })
  return(.out)
}

# Add a layer to the bottom a ggplot object
prepend_layer <- function(x, new_layer) {
  n <- length(x$layers)
  m <- n + 1
  x <- x + new_layer
  x$layers <- c(x$layers[[m]], x$layers[seq(1, n)])
  x
}

dv_pred_ci <- function(
  df,
  x,
  x_lo,
  x_hi,
  y = pm_axis_dv(),
  yname = "value",
  xname = "value",
  xs = list(),
  ys = list(),
  loglog = FALSE,
  scales = c("fixed", "free"),
  ...
) {
  scales <- match.arg(scales)
  if (scales == "fixed") {
    dv_max <- max(c(df[[x_hi]], df[[y]]))
    xs <- c(list(limits = c(0, dv_max)), xs)
    ys <- c(list(limits = c(0, dv_max)), ys)
  }
  p <- dv_pred(df, x, y, yname, xname, xs, ys, loglog, scales, ...)
  p <- prepend_layer(p, geom_errorbarh(
    aes_string(xmin = x_lo, xmax = x_hi), colour = "red", alpha = 0.3
  ))
  return(p)
}
  

