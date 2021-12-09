

# some alternative dispersal kernel functions
dlognormal <- function(x, L, S) (1/((2*pi)^1.5 * S * x^2)) * exp(-(log(x/L))^2 / (2 * S^2))
d2Dt <- function(x, L, S) S / (pi * L * (1 + (x^2 / L))^(S + 1))
dexponential <- function(x, L) dexp(x, rate = L) / (2*pi*x)



#' Construct a neighborhood dispersal probability matrix.
#'
#' @param diameter Neighborhood size, in grid cells (odd integer).
#' @param kernel A dispersal kernel list, e.g. \code{species_template()$kernel}.
#' @param factor Factor by which to increase grid resolution for numerical integration (odd integer).
#' @return A matrix of dispersal probabilities.
#' @export
#' @importFrom stats integrate
#' @importFrom rlang invoke
neighborhood <- function(diameter, kernel, factor = 51){

  # internal functions
  disagg <- function(x, factor) as.matrix(disaggregate(raster(x), factor))
  agg <- function(x, factor) as.matrix(aggregate(raster(x), factor,
                                                 sum, na.rm = T))
  kdf <- function(x){
    kernel$params$x <- x
    invoke(kernel$fun, kernel$params)
  }

  # high-res distance surface
  md <- disagg(matrix(NA, diameter, diameter), factor)
  radius <- (diameter * factor - 1) / 2
  dists <- expand.grid(x = -radius:radius, y = -radius:radius)
  dists <- sqrt(dists$x^2 + dists$y^2) / factor
  dists <- matrix(dists, diameter * factor, diameter * factor)
  dists[dists > radius / factor] <- Inf

  # proportion of propagules falling across neighborhood
  p <- kdf(dists)
  p <- agg(p, factor)
  p <- p / sum(p)

  return(p)
}
