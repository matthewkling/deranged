

# some alternative dispersal kernel functions
dlognormal <- function(x, L, S) (1/((2*pi)^1.5 * S * x^2)) * exp(-(log(x/L))^2 / (2 * S^2))
d2Dt <- function(x, L, S) S / (pi * L * (1 + (x^2 / L))^(S + 1))
dexponential <- function(x, L) dexp(x, rate = L) / (2*pi*x)



# function to aggregate matrices, for internal use
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
#' @param method Either "uniform" (default) or "mixed".
#' @param res Resolution for numerical integration (odd integer).
#' @return A matrix of dispersal probabilities.
#' @export
#' @importFrom stats integrate
#' @importFrom rlang invoke
neighborhood <- function(kernel, diameter = 7, method = "uniform", res = 11){

  kdf <- function(x){
    kernel$params$x <- x
    invoke(kernel$fun, kernel$params)
  }

  # calculate distances and probabilities
  radius <- (diameter * res - 1) / 2
  r <- ifelse(method == "uniform", (res - 1) / 2, 0)
  d <- expand.grid(x = -radius:radius,
                   y = -radius:radius,
                   xo = -r:r,
                   yo = -r:r)
  d$cd <- sqrt(d$x^2 + d$y^2) / res # dist from absolute center
  d$d <- sqrt((d$x - d$xo)^2 + (d$y - d$yo)^2) / res # distance from origin
  d$d <- ifelse(d$cd > (radius / res), Inf, d$d)
  d$p <- kdf(d$d)

  # aggregate across origins
  if(method == "uniform"){
    d <- aggregate(d$p, by = list(xx = d$x, yy = d$y),
                   FUN = mean, na.rm = T)
    d$p <- d$x
  }

  # aggregate across destinations
  p <- agg(matrix(d$p, diameter * res, diameter * res),
           res,
           sum, na.rm = T)

  # normalize
  p <- p / sum(p)
  return(p)
}
