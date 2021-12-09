

# some alternative dispersal kernel functions
dlognormal <- function(x, L, S) (1/((2*pi)^1.5 * S * x^2)) * exp(-(log(x/L))^2 / (2 * S^2))
d2Dt <- function(x, L, S) S / (pi * L * (1 + (x^2 / L))^(S + 1))
dexponential <- function(x, L) dexp(x, rate = L) / (2*pi*x)



# functions to aggregate and disaggregate matrices, for internal use
disagg <- function(x, res){
  y <- matrix(NA, nrow(x) * res, ncol(x) * res)
  xi <- rep(1:nrow(x), each = res)
  xj <- rep(1:ncol(x), each = res)
  for(i in 1:nrow(y)){
    for(j in 1:ncol(y)){
      y[i,j] <- x[xi[i], xj[j]]
    }
  }
  return(y)
}
agg <- function(x, res, fun = sum, ...){
  y <- matrix(NA, nrow(x) / res, ncol(x) / res)
  yi <- rep(1:nrow(y), each = res)
  yj <- rep(1:ncol(y), each = res)
  for(i in 1:nrow(y)){
    for(j in 1:ncol(y)){
      y[i,j] <- fun(x[which(yi == i), which(yj == j)], ...)
    }
  }
  return(y)
}


#' Neighborhood dispersal probability matrix.
#'
#' This function uses numerical integration of a dispersal kernel to construct
#' a matrix representing the probability that dispersal originating in a given
#' grid cell will arrive at each cell across its neighborhood.
#'
#' @param kernel A dispersal kernel list, e.g. \code{species_template()$kernel}.
#' @param diameter Neighborhood size, in grid cells (odd integer).
#' @param origin Either "uniform" (default) or "centroid".
#' @param res Resolution for numerical integration (odd integer).
#' @return A matrix of dispersal probabilities.
#' @export
#' @importFrom stats integrate
#' @importFrom rlang invoke
neighborhood <- function(kernel, diameter = 7, origin = "uniform", res = 11){

  kdf <- function(x){
    kernel$params$x <- x
    invoke(kernel$fun, kernel$params)
  }

  md <- disagg(matrix(NA, diameter, diameter), res)
  radius <- (diameter * res - 1) / 2
  r <- (res - 1) / 2

  if(origin == "centroid"){
    dists <- expand.grid(x = -radius:radius,
                         y = -radius:radius)
    dists$d <- sqrt(dists$x^2 + dists$y^2) / res
    dists$d <- ifelse(dists$d > (radius / res), Inf, dists$d)
    dists$p <- kdf(dists$d)
  }

  if(origin == "uniform"){
    dists <- expand.grid(x = -radius:radius,
                         y = -radius:radius,
                         xo = -r:r,
                         yo = -r:r)
    dists$cd <- sqrt(dists$x^2 + dists$y^2) / res
    dists$d <- sqrt((dists$x - dists$xo)^2 + (dists$y - dists$yo)^2) / res
    dists$d <- ifelse(dists$cd > (radius / res), Inf, dists$d)
    dists$p <- kdf(dists$d)

    dists <- aggregate(dists$p, by = list(xx = dists$x, yy = dists$y),
                       FUN = mean, na.rm = T)
    dists$p <- dists$x
  }

  p <- agg(matrix(dists$p, diameter * res, diameter * res),
           res,
           sum, na.rm = T)
  p <- p / sum(p)
  return(p)
}
