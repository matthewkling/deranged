#include <RcppArmadillo.h>
using namespace Rcpp;


int rbinom_disp(int n,
                double p,
                std::default_random_engine gen) {
  std::binomial_distribution<> d(n, p);
  return d(gen);
}


arma::imat rmultinom_disp(int seeds,
                          arma::mat probs,
                          arma::uvec ind,
                          std::default_random_engine gen) {

  arma::imat y(size(probs), arma::fill::zeros); // post-dispersal counts
  double p = 1; // unallocated probability
  int u = seeds; // unallocated seeds
  arma::uword ii = 0;

  for(arma::uword i = 0; i < probs.n_elem; ++i) {
    ii = ind(i);
    p = accu(probs);
    y(ii) = std::min(rbinom_disp(u, probs(ii) / p, gen), u);
    u = u - y(ii);
    if (u == 0) {
      break;
    }
    probs(ii) = 0;
  }

  return(y);
}


//' Simulate dispersal across a spatial grid
//'
//' @param S A matrix of seed counts across a spatial grid.
//' @param N A neighbor matrix, e.g. produced by \code{neighborhood()}.
//' @param reflect Should dispersers exit the domain (\code{FALSE}) or bounce off the domain boundary (\code{TRUE}, default)?
//' @param rand Randomize dispersal? (default = \code{TRUE})
//' @param seed Integer to seed random number generator.
//' @return A matrix of post-dispersal seed counts of the same dimension as \code{S}.
//' @export
// [[Rcpp::export]]
arma::mat disperse(arma::mat S,
                   arma::mat N,
                   bool reflect = true,
                   bool rand = true,
                   int seed = 1) {

  std::default_random_engine gen(seed); // initialize random number generator

  arma::uvec Ni = arma::sort_index(N, "descent"); // order to evaluate neighbors

  int r = (N.n_rows - 1) / 2; // window radius
  arma::mat T(S.n_rows + r * 2, S.n_cols + r * 2, arma::fill::zeros); // padded grid

  for(arma::uword a = 0; a < S.n_rows; ++a) {
    for(arma::uword b = 0; b < S.n_cols; ++b) {

      if (rand) {
        T.submat(a, b, a + r * 2, b + r * 2) =
          T.submat(a, b, a + r * 2, b + r * 2) +
          rmultinom_disp(S(a, b), N, Ni, gen);
      } else {
        T.submat(a, b, a + r * 2, b + r * 2) =
          T.submat(a, b, a + r * 2, b + r * 2) +
          S(a, b) * N;
      }

    }
  }

  if (reflect) {
    for(int i = 0; i < r; ++i){
      T.row(r * 2 - 1 - i) = T.row(r * 2 - 1 - i) + T.row(i);
      T.col(r * 2 - 1 - i) = T.col(r * 2 - 1 - i) + T.col(i);
      T.row(T.n_rows - (r * 2 - 1 - i) - 1) =
        T.row(T.n_rows - (r * 2 - 1 - i) - 1) + T.row(T.n_rows - 1 - i);
      T.col(T.n_cols - (r * 2 - 1 - i) - 1) =
        T.col(T.n_cols - (r * 2 - 1 - i) - 1) + T.col(T.n_cols - 1 - i);
    }
  }

  return T.submat(r, r, S.n_rows + r - 1, S.n_cols + r - 1);
}

