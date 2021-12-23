

#' Generate data structure for species-level parameters.
#'
#' @param n_env Number of environmental variables (integer).
#' @param names Demographic stage names (character vector).
#' @return A list of demographic and dispersal parameter objects.
#' @export
species_template <- function(n_env = 1,
                             names = c("s", "j", "a")){

  list(alpha = matrix(0, 3, 3,
                      dimnames = list(names, names)),
       beta = array(0, c(3, 3, 3),
                    dimnames = list(names, names, names)),
       gamma = array(0, c(3, 3, n_env),
                     dimnames = list(names, names, paste0("v", 1:n_env))),
       fecundity = setNames(rep(0, length(names)), names),
       kernel = list(fun = dlognormal,
                     params = list(L = 1, S = 1)))
}


#' Generate example data structure for spatial data.
#'
#' @param n_row Number of rows in spatial grid.
#' @param n_col Number of columns in spatial grid.
#' @param n_env Number of environmental variables.
#' @param n_steps Number of time steps for environmental variables.
#' @param names Demographic stage names (character vector).
#' @return A list of demographic and dispersal parameter objects.
#' @export
landscape_template <- function(n_row = 10, n_col = 10,
                               n_steps = 1, n_env = 1,
                               names = c("s", "j", "a")){

  # environmental data -- 4d array (space, space, variable, time)
  e <- array(0, c(n_row, n_col, n_env, n_steps),
             dimnames = list(NULL, NULL, paste0("v", 1:n_env), paste0("t", 1:n_steps)))

  # initial population rasters -- 3d array (space, space, stage)
  n <- array(0, c(n_row, n_col, length(names)),
             dimnames = list(NULL, NULL, names))

  list(e = e,
       n = n)
}


#' Run a range simulation
#'
#' @param sp Species parameter list, following \code{species_template()}.
#' @param ls Landscape spatial data list, following \code{landscape_template()}.
#' @param n_steps Number of time steps to simulate (integer).
#' @param randomize Should demography and dispersal be randomized (logical)?
#' @param reflect Should dispersers bounce off domain boundary (logical)?
#' @param record Index of age class to record and return (integer).
#' @param ... Further arguments passed to \code{neighborhood()}.
#' @return An array of population values over space and time, for the class specified in \code{record}.
#' @export
#' @importFrom utils setTxtProgressBar txtProgressBar
simulate <- function(sp,
                     ls,
                     n_steps = 100,
                     randomize = TRUE,
                     reflect = TRUE,
                     record = 3,
                     seed = 1,
                     ...){

  sim(N = ls$n,
      env = lapply(1:dim(ls$e)[4], function(i) array(ls$e[,,,i], dim(ls$e)[1:3])),
      alpha = sp$alpha, beta = sp$beta, gamma = sp$gamma, fecundity = sp$fecundity,
      nb = neighborhood(sp$kernel, ...),
      nsteps = n_steps,
      rand = randomize,
      reflect = reflect,
      record = record - 1,
      seed = seed)
}

