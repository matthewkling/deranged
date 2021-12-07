

#' Generate data structure for species-level parameters.
#'
#' @param n_env Number of environmental variables (integer).
#' @param names Class names (character vector).
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
#' @param n_steps Number of time steps for environmental variables.
#' @param n_env Number of environmental variables.
#' @param names Class names (character vector).
#' @return A list of demographic and dispersal parameter objects.
#' @export
landscape_template <- function(n_row = 10, n_col = 10,
                               n_steps = 1, n_env = 1,
                               names = c("s", "j", "a")){

  # environmental data -- 4d array (space, space, time, variable)
  e <- array(0, c(n_row, n_col, n_steps, n_env),
             dimnames = list(NULL, NULL, paste0("t", 1:n_steps), paste0("v", 1:n_env)))

  # initial population rasters -- 3d array (space, space, time)
  n <- array(0, c(n_row, n_col, 3),
             dimnames = list(NULL, NULL, names))

  list(e = e,
       n = n)
}





#' Run a multi-year range simulation
#'
#' @param sp Species parameter list, following \code{species_template()}.
#' @param ls Landscape spatial data list, following \code{landscape_template()}.
#' @param diameter Neighborhood size (integer).
#' @param n_steps Number of time steps to simulate (integer).
#' @param randomize Should demography and dispersal be randomized (logical)?
#' @param reflect Should dispersers bounce off domain boundary (logical)?
#' @param record Index of age class to record and return (integer).
#' @return An array of population values over space and time, for the class specified in \code{record}.
#' @export
#' @importFrom utils setTxtProgressBar txtProgressBar
simulate <- function(sp,
                     ls,
                     diameter = 7,
                     n_steps = 100,
                     randomize = TRUE,
                     reflect = TRUE,
                     record = 3){

  neighbors <- neighborhood(diameter, sp$kernel)

  n <- as.array(ls$n)

  d <- array(NA, c(dim(n)[1:2], n_steps + 1),
             dimnames = list(1:dim(n)[1],
                             1:dim(n)[2],
                             0:n_steps))
  d[ , , 1] <- n[ , , record]

  e <- as.array(ls$e)
  if(dim(e)[3] == n_steps){ ei <- 1:n_steps}else{ei <- rep(1L, n_steps)}

  pb <- txtProgressBar(min = 0, max = n_steps, initial = 0, style = 3)
  for(i in 1:n_steps){
    d[ , , i + 1] <- n[ , , record]
    n <- transition(n,
                    E = array(e[, , ei[i], ], dim(e)[c(1, 2, 4)]),
                    alpha = sp$alpha, beta = sp$beta, gamma = sp$gamma,
                    rand = randomize,
                    seed = sample(1e8, 1))
    n[, , 1] <- n[, , 1] + disperse(reproduce(n, f = sp$fecundity),
                                    neighbors,
                                    reflect = reflect,
                                    rand = randomize,
                                    seed = sample(1e8, 1))
    setTxtProgressBar(pb, i)
  }

  return(d)
}

