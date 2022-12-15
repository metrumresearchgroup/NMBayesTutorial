### Helper functions for modeling with bbr ####  ----------------------------------

#' Helper to open control stream for editing in Rstudio
#' @param .mod a `bb{.model_type}_model` object
edit_model <- function(.mod){
  .mod %>%
    get_model_path() %>% 
    file.edit()
}

#' Check tags against external glossary
#' 
#' Checks whether any of `tags` are _not_ contained in
#' the list or vector passed to `glossary`
audit_tags <- function(tags, glossary) {
  UseMethod("audit_tags")
}

#' @describeIn audit_tags Takes a character vector of tags and 
#' checks if any are not contained in `glossary`
audit_tags.character <- function(tags, glossary = yaml::read_yaml("tags.yaml")) {
  tags_bool <- tags %in% glossary
  
  if (any(!tags_bool)) {
    warning(paste(
      "The following are not valid tags:",
      paste(tags[which(!tags_bool)], collapse = ", ")
    ), call. = FALSE)
    return(invisible(FALSE))
  }
  return(invisible(TRUE))
}

#' @describeIn audit_tags Takes a tibble of class `bbi_run_log_df` (the output of [bbr::run_log()])
#' and checks if any of the models contain tags that are not contained in `glossary`
audit_tags.bbi_run_log_df <- function(tags, glossary = yaml::read_yaml("tags.yaml")) {
  tags_bool <- purrr::map(log_df$tags, function(.t) {
    .t %in% glossary
  } )
  
  bad_mods <- purrr::map_lgl(tags_bool, ~any(!.x))
  
  if (any(bad_mods)) {
    bad_msg <- purrr::map_chr(which(bad_mods), function(.m) {
      mod_path <- log_df$absolute_model_path[.m]
      bad_tags <- log_df$tags[[.m]][!tags_bool[[.m]]]
      paste0("  ", mod_path, ": ", paste(bad_tags, collapse = ", "))
    })
    
    warning(paste(c(
      "The following models have invalid tags:",
      bad_msg
    ), collapse = "\n"))
  }
  
  invisible(!bad_mods)
}


#' Compare objective function values between models
#'
#' Takes a bbi_summary_log_df and character vector of models 
#' and filters to only those models. 
#' Returns tibble with run, ofv, tags,
#' and any columns passed to ...
compare_ofv <- function(models, log_df, ...) {
  log_df %>%
    filter(run %in% models) %>%
    select(run, ofv, tags, ...) %>%
    collapse_to_string(tags, ...) %>%
    knitr::kable()
}

#' helper to check if model or data is out of date
#' 
#' Checks the `model_has_changed` and `data_has_changed` 
#' in a `bbi_config_log_df` and returns the offending rows
#' if there are any.
check_up_to_date <- function(log_df, .kable = TRUE) {
  
  if (!inherits(log_df, "bbi_config_log_df")) {
    log_df <- log_df %>% add_config()
  }
  
  # get md5 digests and compare to those from config_log() columns
  problem_df <- log_df %>%
    filter(data_has_changed | model_has_changed)
  
  if (nrow(problem_df) != 0) {
    cat("## Warning: data or model has changed\nIn the following models, either the data or model have changed since the time they were run.\n")
    
    problem_df <- problem_df %>%
      select(run, data_path, data_has_changed, model_has_changed) 
    
    if (isTRUE(.kable)) {
      problem_df <- knitr::kable(problem_df)
    }
    
   problem_df
  } else {
    cat("All models and data are up to date.\n")
  }
}


#' Update control stream run numbers
#' 
#' Helper to update the run number in the new control stream of a model
#' created with [bbr::copy_model_from()]. See details section.
#' 
#' @details This function updates all occurances of the run number from
#' the parent model and replaces them with run number from the 
#' new model. 
#' 
#' The function relies on the assumption that `.mod` was created by 
#' `copy_model_from()` and therefore will have the parent model number as
#' the first entry in `.mod$based_on`.
#' 
#' It will look for that parent model number before all strings 
#' passed to `.suffixes` (e.g. `{parent_mod}.MSF`, etc.)
#' replace it with the `get_model_id(.mod)` wherever found.
#' `.suffixes` defaults to the following:
#' 
#' * `.MSF`
#' * `.EXT`
#' * `.CHN`
#' * `.tab`
#' * `par.tab`
#' 
#' All matches are _not_ case sensitive, meaning `{parent_mod}.MSF` and 
#' `{parent_mod}.msf` would both be replaced by `{get_model_id(.mod)}.MSF`.
#' 
#' @param .mod The `bbi_nonmem_model` object associated with the control
#'   stream that will be modified.
#' @param .suffixes Character vector of suffixes to be matched for replacement.
#'   Matching is case insensitive.
update_run_number <- function(
  .mod,
  .suffixes = c(
    '.MSF',
    '.EXT',
    '.CHN',
    '.tab',
    'par.tab'
  )
){
  runno <- get_model_id(.mod)
  modelfile <- get_model_path(.mod)
  based_on <- .mod$based_on[1]
  message(glue::glue("replacing {based_on} with {runno} in {modelfile}"))
  
  ## edit text of new model file
  txt <- suppressSpecificWarning(readLines(modelfile), .regexpr = "incomplete final line")
  for (.s in .suffixes) {
    txt <- gsub(
      paste0(based_on, .s), 
      paste0(runno, .s), 
      txt, 
      ignore.case = T
    )
  }
  
  ## write updated model file
  cat(txt, file=modelfile, sep='\n')
  
  ## return model to make this a pipeable function
  return(.mod)
}


#' Helper to show the diff between two control streams
#' 
#' @param model_A `bbi_{.model_type}_model` object to compare
#' @param model_B `bbi_{.model_type}_model` object to compare
#' @param render_in_md If `FALSE`, the default, will render the diff
#'   in the viewer panel. This is for interactive viewing.
#'   To use this function to render a diff in Rmd, set this to `TRUE`
#'   and set the relevant block to `results = "asis"`.
model_diff <- function(model_A, model_B, render_in_rmd = FALSE) {
  model_A <- get_model_path(model_A)
  model_B <- get_model_path(model_B)
  
  suppressSpecificWarning({
    if (isTRUE(render_in_rmd)) {
      diffobj::diffFile(
        model_A,
        model_B,
        mode = "sidebyside",
        format = "html",
        style = list(html.output = "diff.w.style"),
        interactive = F
      )
    } else {
      diffobj::diffFile(
        model_A,
        model_B,
        mode = "sidebyside"
      )
    }
  }, .regexpr = "incomplete final line") 
}

#' Helper to show the diff between tags
#' 
#' @param model_A `bbi_{.model_type}_model` object to compare
#' @param model_B `bbi_{.model_type}_model` object to compare
#' @param .print If `TRUE`, the default, will print to console. 
#' @return A named list with an element for both models containing
#'   the tags that are in that model but _not_ the other model. 
#'   Returns invisibly since `.print` is `TRUE` by default.
tags_diff <- function(model_A, model_B, .print = TRUE) {
  diff_A <- setdiff(model_A$tags, model_B$tags)
  diff_B <- setdiff(model_B$tags, model_A$tags)
  
  id_A <- get_model_id(model_A)
  id_B <- get_model_id(model_B)
  
  if (isTRUE(.print)) {
    cat(paste(
      glue::glue("In {id_A} but not {id_B}:\t{paste(diff_A, collapse = ', ')}\n"),
      glue::glue("In {id_B} but not {id_A}:\t{paste(diff_B, collapse = ', ')}\n"),
      sep = "\n"
    ))
  }
  
  .out <- list()
  .out[[id_A]] <- diff_A
  .out[[id_B]] <- diff_B
  
  invisible(.out)
}

#' Run Bayes chains
#' 
#' Run multiple chains of a Bayes model after initial estimates have been
#' generated
#'
#' @param .model_dir path to directory containing model
#' @param .run run name
#' @param .bbi_args list of arguments to be passed to `submit_model()`
run_chains <- function(.model_dir, .run, .mode = "sge", .bbi_args) {
  mod <- read_model(file.path(.model_dir, .run))
  ctl <- read_lines(get_model_path(mod))
  
  row_bayes <- str_detect(ctl, "METHOD=BAYES|METHOD=NUTS")
  est_bayes <- ctl[row_bayes]
  est_bayes <- str_replace(est_bayes, "^;", "")
  ctl[row_bayes] <- est_bayes
  
  row_table <- str_detect(ctl, ";\\s*\\$TABLE")
  block_table <- ctl[row_table]
  block_table <- str_replace(block_table, "^;", "")
  ctl[row_table] <- block_table
  
  row_chain <- str_detect(ctl, "METHOD=CHAIN")
  est_chain <- ctl[row_chain]
  n_chain <- as.numeric(str_extract(est_chain, "(?<=NSAMPLE=)[0-9]+"))
  est_chain <- str_replace(est_chain, "NSAMPLE=[0-9]+", "NSAMPLE=0")
  est_chain <- str_replace(est_chain, "FILE=", "FILE=../")
  
  row_data <- str_detect(ctl, "\\$DATA")
  data_record <- ctl[row_data]
  ctl[row_data] <- str_replace(data_record, "\\$DATA\\s+", "$DATA ../")
  
  row_extrasend <- str_detect(ctl, "extrasend")
  ctl[row_extrasend] <- str_replace(ctl[row_extrasend], "extrasend", "../extrasend")
  
  walk(seq_len(n_chain), function(.chain) {
    #cat(.chain, "\n")
    est_chain_i <- str_replace(est_chain, "ISAMPLE=0", glue::glue("ISAMPLE={.chain}"))
    #est_chain_i <- str_replace(est_chain_i, "SEED=[0-9]+", glue::glue("SEED={.chain}"))
    est_bayes_i <- str_replace(est_bayes, "SEED=[0-9]+", glue::glue("SEED={.chain}"))
    #cat(est_chain_i, "\n")
    ctl_i <- ctl
    ctl_i[row_chain] <- est_chain_i
    ctl_i[row_bayes] <- est_bayes_i
    write_lines(ctl_i, file.path(
      .model_dir,
      glue::glue("{.run}/{.run}_{.chain}.ctl"))
    )
    
    mod <- new_model(
      #glue::glue("./{.run}/{.run}.{.chain}.yaml"),
      file.path(.model_dir, .run, glue::glue("{.run}_{.chain}")),
      .description = glue::glue("Chain {.chain}"),
      .overwrite = TRUE
    )
    
    proc <- submit_model(
      mod,
      #.directory = file.path(.model_dir, .run),
      .bbi_args = .bbi_args,
      .mode = .mode,
      .config_path = file.path(.model_dir, "bbi.yaml"),
      .wait=FALSE
      #.dry_run = FALSE
    )
  })
}

#' Return a vector of absolute paths to table files to look for
#' 
#' @param mod the `bbi_nonmem_model` object
#' 
#' @details
#' Called by redataset
#' 
tabfiles <- function(mod) {
  output_dir <- get_output_dir(mod)
  c(
    file.path(output_dir, "PAR_TAB"), 
    file.path(output_dir, "TAB"),
    build_path_from_model(mod, "par.tab"), 
    build_path_from_model(mod, ".tab")
    
  )
}


#' Return a single data frame with model output and input data
#' 
#' @param mod the `bbi_nonmem_model` object or a path to a NONMEM run
#' @param join character column name to use to join tab files
#' @param files absolute paths to table files to try to join
#' @param more absolute paths to other files to try to join
#' @param verbose if TRUE, messages are printed as files are read
#' 
#' @details
#' Required packges: data.table, dplyr, bbr. 
#' 
#' Note that you can get the absolute path to the model output directory 
#' with `bbr::get_output_dir(mod)` if needed to build paths for passing to 
#' `files` or `more`.
#' 
redataset <- function(mod, join = "NUM", files = tabfiles(mod), more = NULL, 
                      verbose = TRUE) {
  stopifnot(requireNamespace("data.table"))
  stopifnot(requireNamespace("dplyr"))
  stopifnot(requireNamespace("bbr"))
  if (!inherits(mod, "bbi_nonmem_model")) {
    mod <- read_model(mod)
  }
  files <- c(files,more)
  verbose <- isTRUE(verbose)
  sum <- bbr::model_summary(mod)
  nid <-  sum$run_details$number_of_patients
  nrec <- sum$run_details$number_of_data_records
  data_loc <- get_data_path(sum)
  if(verbose) message("data file: ", basename(data_loc))

  data <- data.table::fread(data_loc, na = '.', data.table = FALSE)
  data <- dplyr::as_tibble(data)
  if(verbose) {
    message("  rows: ", nrow(data))
    message("  cols: ", ncol(data))
  }
  chk <- file.exists(files)
  files <- files[chk]
  if(length(files)==0) {
    message("  zero table files found; returning")
    return(data)
  }
  if(!all(join %in% names(data))) {
    stop("couldn't find '", join, "' column in 'data'")  
  }
  leading_cols <- names(data)
  for(p in files) {
    tab <- data.table::fread(p, skip = 1, data.table = FALSE, na ='.')
    tab <- dplyr::as_tibble(tab)
    has_id <- "ID" %in% names(tab)
    if(!(nrow(tab) %in% c(nrec,nid))) {
      message("\ntable file: ", basename(p), " (skipped)")
      message("  rows: ", nrow(tab))
      message("  hasid: ", has_id)
      message("  tabids: ", length(unique(tab[["ID"]])))
      message("  runids: ", nid)
      next;
    }
    if(nrow(tab) == nid && has_id) {
      if(verbose) {
        message("\nfirstonly file: ", p)
        message("  rows: ", nrow(tab))
        message("  cols: ", ncol(tab))
      }
      keep <- c("ID",setdiff(names(tab), names(data)))
      if((ncol(tab) - length(keep) - 1) > 0) {
        if(verbose) {
          message(
            "  droppping: ", 
            ncol(tab) - length(keep) - 1, 
            " columns in common"
          )
        }
      }
      data <- dplyr::left_join(data, tab[,keep], by = "ID", suffix = c("", "y"))
      next;
    } # end ID join
    if(nrow(tab) != nrec) {
      message("\ntable file: ", basename(p), " (skipped)")
      message("  rows: ", nrow(tab))
      message("  hasid: ", has_id)
      message("  tabids: ", length(unique(tab[["ID"]])))
      message("  runids: ", nid)
      next;
    }
    if(verbose) message("\ntable file: ", basename(p))
    new_cols <- setdiff(names(tab), names(data))
    keep <- c(join,new_cols)
    drop <- setdiff(names(tab), new_cols)
    drop <- drop[!(drop %in% join)]
    if(verbose) {
      message("  rows: ", nrow(tab))
      message("  cols: ", length(new_cols), " new")
    }
    for(d in drop) {
      message("  drop: ", d,  " common")  
    }
    data <- dplyr::left_join(tab[,keep], data, by = join)
  }  
  if(verbose) {
    message("\nfinal stats:")
    message("  rows: ", nrow(data))
    message("  cols: ", ncol(data))
  }
  dplyr::select(data, dplyr::all_of(leading_cols), dplyr::everything())
}
