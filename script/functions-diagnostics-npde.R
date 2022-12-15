##' Purpose: 
##' define functions that are called repeatedly while producing diagnostics
##' for model runs
##' 


##' map_wrap_eta_cont
##' Purpose: plot ETA vs all continuous covariates
##' @param  .map_etas character: name of ETA
##' @param  .co       string: continuous covariate list
##' @param  .id       dataframe: includes one record per id
##' @param  .ncol     numeric: number of columns in plot
map_wrap_eta_cont <- function(.map_etas,.co,.id, 
                              .ncol = 2) {
  .p <- wrap_eta_cont(
    .id,
    y = .map_etas,
    x = .co,
    use_labels = TRUE,
    ncol = .ncol, scales= "free_x"
  )
}
##' map_eta_cat
##' Purpose: plot all ETAs vs a categorical covariate
##' @param  .map_ca   character: name of a categorical covariate
##' @param  .etas     string: ETA list
##' @param  .id       dataframe: includes one record per id
map_eta_cat <- function(.map_ca, .etas, .id, .rot = 45) {
  .p <- eta_cat(.id, x = .map_ca, y = .etas) %>% 
    ## CHECK: depending on the labels, this may need to be changed 
    purrr::map(~.x+rot_x(.rot)) %>% 
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
      ),
      skip=1   ## Included this for PD models
    ) %>% 
      as_tibble() %>% 
      mutate(chain = .chain)
  })
}

decorr.chol <- function(x) {
  xerr<-0
  xmat<-try(chol(posmat(x)))
  if(is.numeric(xmat)) {
    ymat<-try(solve(xmat))
    if(!is.numeric(ymat))
      xerr<-2
  } else
    xerr<-1
  return(list(y=ymat,xerr=xerr))
}

decorr.inverse <- function(x) {
  xerr<-0
  var.eig<-eigen(posmat(x))
  xmat<-try(var.eig$vectors %*% diag(sqrt(var.eig$values)) %*% solve(var.eig$vectors))
  if(is.numeric(xmat)) {
    ymat<-try(solve(xmat))
    if(!is.numeric(ymat))
      xerr<-2
  } else
    xerr<-1
  return(list(y=ymat,xerr=xerr))
}

decorr.polar <- function(x) {
  xerr<-0
  xmat<-try(chol(posmat(x)))
  if(is.numeric(xmat)) {
    svdec<-svd(xmat)
    umat<-svdec$u %*% t(svdec$v)
    vmat<-t(umat) %*% xmat
    ymat<-try(solve(vmat))
    if(!is.numeric(ymat))
      xerr<-2
  } else
    xerr<-1
  return(list(y=ymat,xerr=xerr))
}

is.square <- function(x,...) dim(x)[[1]] == dim(x)[[2]]

posmat <- function(x,...) {
  if(any(diag(x) <=0)) stop("matrix cannot be made positive-definite")
  if(!is.square(x))stop('x is not square')
  sign <- sign(x)
  x <- abs(x)
  characteristic <- trunc(log(x,10))
  mantissa <- log(x,10) - characteristic
  scale <- 10^characteristic
  digits <- 10^mantissa * 1e5
  diagonal <- round(diag(digits))
  digits <- floor(digits)
  diag(digits) <- diagonal
  digits <- digits/1e5
  x <- sign * scale * digits
  diagonal <- diag(x)
  y <- 0.97 * x
  diag(y) <- diagonal
  if(det(x)>0) x else posmat(y)
}


assignInNamespace('decorr.chol', decorr.chol, "npde")
assignInNamespace('decorr.inverse', decorr.inverse, "npde")
assignInNamespace('decorr.polar', decorr.polar, "npde")

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
