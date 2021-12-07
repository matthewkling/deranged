

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
