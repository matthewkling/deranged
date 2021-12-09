

# some alternative dispersal kernel functions
dlognormal <- function(x, L, S) (1/((2*pi)^1.5 * S * x^2)) * exp(-(log(x/L))^2 / (2 * S^2))
d2Dt <- function(x, L, S) S / (pi * L * (1 + (x^2 / L))^(S + 1))
dexponential <- function(x, L) dexp(x, rate = L) / (2*pi*x)



#' Neighborhood dispersal probability matrix.
#'
#' This function uses numerical integration of a dispersal kernel to construct
#' a matrix representing the probability that dispersal originating in a given
#' grid cell will arrive at each cell across its neighborhood.
#'
#' @param diameter Neighborhood size, in grid cells (odd integer).
#' @param kernel A dispersal kernel list, e.g. \code{species_template()$kernel}.
#' @param origin Either "uniform" (default) or "centroid".
#' @param res Resolution for numerical integration (odd integer).
#' @return A matrix of dispersal probabilities.
#' @export
#' @importFrom stats integrate
#' @importFrom rlang invoke
neighborhood <- function(diameter, kernel, origin = "uniform", res = 11){

  disagg <- function(x, res) as.matrix(disaggregate(raster(x), res))

  agg <- function(x, res) as.matrix(aggregate(raster(x), res, sum, na.rm = T))

  kdf <- function(x){
    kernel$params$x <- x
    invoke(kernel$fun, kernel$params)
  }

  md <- disagg(matrix(NA, diameter, diameter), res)
  radius <- (diameter * res - 1) / 2
  r <- (res - 1) / 2

  if(origin == "centroid"){
    dists <- expand_grid(x = -radius:radius,
                         y = -radius:radius) %>%
      mutate(d = sqrt(x^2 + y^2) / res,
             d = ifelse(d > (radius / res), Inf, d),
             p = kdf(d))
  }

  if(origin == "uniform"){
    dists <- expand_grid(x = -radius:radius,
                         y = -radius:radius,
                         xo = -r:r,
                         yo = -r:r) %>%
      mutate(cd = sqrt(x^2 + y^2) / res,
             d = sqrt((x-xo)^2 + (y-yo)^2) / res,
             d = ifelse(cd > (radius / res), Inf, d),
             p = kdf(d)) %>%
      group_by(x, y) %>%
      summarize(p = mean(p, na.rm = T)) %>%
      ungroup()
  }

  p <- dists$p %>%
    matrix(diameter * res, diameter * res) %>%
    agg(res)
  p <- p / sum(p)
  return(p)
}
