#include <RcppArmadillo.h>
using namespace Rcpp;


double rbinom_mms(arma::mat n,
                  arma::mat p,
                  int seed = 1) {
  std::default_random_engine gen(seed);
  arma::vec y(n.n_elem);
  for(arma::uword i = 0; i < n.n_elem; ++i) {
    std::binomial_distribution<> d(n(i), p(i));
    y(i) = d(gen);
  }
  return sum(y);
}

arma::mat rbinom_mmm(arma::mat n,
                     arma::mat p,
                     int seed = 1) {
  std::default_random_engine gen(seed);
  arma::mat y(n.n_rows, n.n_cols);
  for(arma::uword i = 0; i < n.n_elem; ++i) {
    std::binomial_distribution<> d(n(i), p(i));
    y(i) = d(gen);
  }
  return y;
}

//' Seed dispersal across a spatial grid
//'
//' @param S A matrix of seed counts across a spatial grid.
//' @param N A neighbor matrix, e.g. produced by \code{neighborhood()}.
//' @param rand Randomize dispersal? (Boolean, default = TRUE).
//' @param seed Integer to seed random number generator.
//' @return A matrix of post-dispersal seed counts of the same dimension as \code{S}.
//' @export
// [[Rcpp::export]]
arma::mat disperse(arma::mat S,
                   arma::mat N,
                   bool rand = true,
                   int seed = 1) {
  int r = (N.n_rows - 1) / 2;
  arma::mat T(S.n_rows + r * 2, S.n_cols + r * 2, arma::fill::zeros);
  T.submat(r, r, S.n_rows + r - 1, S.n_cols + r - 1) = S;
  arma::mat U(S.n_rows, S.n_cols, arma::fill::zeros);

  for(arma::uword a = 0; a < S.n_rows; ++a) {
    for(arma::uword b = 0; b < S.n_cols; ++b) {

      if (rand) {
        U(a, b) = rbinom_mms(T.submat(a, b, a + r * 2, b + r * 2), N, seed);
        ++seed;
      } else {
        U(a, b) = accu(T.submat(a, b, a + r * 2, b + r * 2) % N);
      }

    }
  }
  return U;
}


//' Perform a randomized stage-based demographic transition
//'
//' @param N A 3-D array of population numbers for each life stage, over a spatial grid (x, y, class).
//' @param E A 4-D array of environmental data (x, y, time, variable).
//' @param alpha A matrix of transition intercepts (from, to).
//' @param beta A 3-D array of density dependence effects (from, to, modifier).
//' @param gamma A 3-D array of environmental effects (from, to, variable).
//' @param fecundity Multiplier applied to first life stage after transition.
//' @param rand Randomize transitions instead of using matrix multiplication? (Boolean, default = TRUE).
//' @param seed Integer to seed random number generator.
//' @return A 3-D array of population numbers for each life stage.
//' @export
// [[Rcpp::export]]
arma::cube transition(arma::cube N,
                      arma::cube E,
                      arma::mat alpha,
                      arma::cube beta,
                      arma::cube gamma,
                      double fecundity = 1,
                      bool rand = true,
                      int seed = 1) {

  arma::cube NN(size(N), arma::fill::zeros);

  arma::mat p(N.n_rows, N.n_cols, arma::fill::zeros);
  double m = 0;

  for(arma::uword t = 0; t < alpha.n_cols; ++t) { // target (to)
    for(arma::uword s = 0; s < alpha.n_rows; ++s) { // source (from)

      if (alpha(s, t) +
          accu(beta.tube(s, t)) +
          accu(gamma.tube(s, t)) == 0) {
        continue;
      }

      p.fill(alpha(s, t)); // intercept

      for(arma::uword d = 0; d < N.n_slices; ++d){ // density dependence
        m = beta(s, t, d);
        if (m == 0) {
          continue;
        }
        p = p + N.slice(d) * m;
      }

      for(arma::uword e = 0; e < E.n_slices; ++e){ // environmental effects
        m = gamma(s, t, e);
        if (m == 0) {
          continue;
        }
        p = p + E.slice(e) * m;
      }

      if (rand) {
        NN.slice(t) = NN.slice(t) + rbinom_mmm(N.slice(s), clamp(p, 0, 1), seed);
        ++seed;
      } else {
        NN.slice(t) = NN.slice(t) + N.slice(s) % clamp(p, 0, 1);
      }

    }
  }

  if (fecundity != 1) {
    NN.slice(0) = NN.slice(0) * fecundity;
  }

  return NN;
}

