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




