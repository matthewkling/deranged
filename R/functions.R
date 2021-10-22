

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


# some alternative dispersal kernel functions
dlognormal <- function(x, L, S) (1/((2*pi)^1.5 * S * x^2)) * exp(-(log(x/L))^2 / (2 * S^2))
d2Dt <- function(x, L, S) S / (pi * L * (1 + (x^2 / L))^(S + 1))
dexponential <- function(x, L) dexp(x, rate = L) / (2*pi*x)


#' Construct a neighborhood dispersal probability matrix.
#'
#' @param diameter Neighborhood size, in grid cells (integer).
#' @param kernel A dispersal kernel list, e.g. \code{species_template()$kernel}.
#' @return A matrix of dispersal probabilities.
#' @export
#' @importFrom stats integrate
#' @importFrom rlang invoke
neighborhood <- function(diameter, kernel){

  # neighbor distances (units = original grid cells)
  if(diameter %% 2 != 1) stop("diameter must be an odd integer")
  radius <- (diameter - 1) / 2
  dists <- expand.grid(x = -radius:radius, y = -radius:radius)
  dists <- sqrt(dists$x^2 + dists$y^2)
  dists <- matrix(dists, diameter, diameter)
  dists[dists > radius] <- Inf

  # proportion of propagules falling across neighborhood
  kdf <- function(x){
    kernel$params$x <- x
    invoke(kernel$fun, kernel$params)
  }
  kd <- function(x){
    d <- kdf(x)
    p <- integrate(function(x) kdf(x)*2*pi*x,
                   0, .5)$value # fraction of propagules not leaving cell
    r <- (nrow(x) + 1)/2
    d[r, r] <- 0
    d <- d / sum(d) * (1 - p)
    d[r, r] <- p
    return(d)
  }
  kd(dists)
}


#' Run a multi-year range simulation
#'
#' @param sp Species parameter list, following \code{species_template()}.
#' @param ls Landscape spatial data list, following \code{landscape_template()}.
#' @param diameter Neighborhood size (integer).
#' @param n_steps Number of time steps to simulate (integer).
#' @param record Index of age class to record and return (integer).
#' @param randomize Should demography and dispersal be randomized (logical)?
#' @return An array of population values over space and time, for the class specified in \code{record}.
#' @export
#' @importFrom utils setTxtProgressBar txtProgressBar
simulate <- function(sp,
                     ls,
                     diameter = 7,
                     n_steps = 100,
                     record = 3,
                     randomize = TRUE){

  neighbors <- neighborhood(diameter, sp$kernel)

  n <- as.array(ls$n)

  d <- array(NA, c(dim(n)[1:2], n_steps),
             dimnames = list(1:dim(n)[1],
                             1:dim(n)[2],
                             1:n_steps))

  e <- as.array(ls$e)
  if(dim(e)[3] == n_steps){ ei <- 1:n_steps}else{ei <- rep(1L, n_steps)}

  pb <- txtProgressBar(min = 0, max = n_steps, initial = 0, style = 3)
  for(i in 1:n_steps){
    d[ , , i] <- n[ , , record]
    n <- transition(n,
                    E = array(e[, , ei[i], ], dim(e)[c(1, 2, 4)]),
                    alpha = sp$alpha, beta = sp$beta, gamma = sp$gamma,
                    rand = randomize)
    n[, , 1] <- disperse(n[, , 1], neighbors, rand = randomize)
    setTxtProgressBar(pb, i)
  }

  return(d)
}

