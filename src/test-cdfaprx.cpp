#include "cdfaprx.h"
#include <testthat.h>
#include <limits>
#include "threat-safe-random.h"

context("restrictcdf unit tests") {
  test_that("cdf<likelihood> gives similar output to R") {
    /*
     set.seed(1)
     n <- 4
     mean <- rnorm(n)
     sigma <- drop(rWishart(1L, 2L * n, diag(n)))

     lower <- c(-Inf, -1 ,  1, -Inf)
     upper <- c(   1, Inf,  3,    2)
     mean  <- round(mean , 3)
     sigma <- round(sigma, 3)

     library(mvtnorm)
     prob <- pmvnorm(lower, upper, mean, sigma = sigma,
     algorithm = GenzBretz(abseps = 1e-9, maxpts = 1000000L))

     dput(mean)
     dput(sigma)
     dput(prob)
     */
    std::vector<unsigned> seeds = { 1L };
    parallelrng::set_rng_seeds(seeds);
    constexpr double const Inf = std::numeric_limits<double>::infinity();

    arma::vec lower, upper, mean;
    arma::mat sigma;

    lower << -Inf   << -1   << 1      << -Inf;
    upper << 1      << Inf  << 3      << 2;
    mean  << -0.626 << 0.18 << -0.836 << 1.595;

    sigma << 8.287 << -0.848 << -0.879 << -1.788 << -0.848 << 3.581
          << 2.916 << -3.957 << -0.879 << 2.916 << 7.361 << -0.648
          << -1.788 << -3.957 << -0.648 << 11.735;
    sigma.reshape(4L, 4L);

    double constexpr E_prop(0.0693596863013216);
    constexpr double abs_eps{E_prop * 1e-3};
    pedmod::cdf<pedmod::likelihood>::alloc_mem(4, 1);
    pedmod::likelihood::alloc_mem(4, 1, 8);
    pedmod::likelihood func;

    {
      auto res = pedmod::cdf<pedmod::likelihood>(
        func, lower, upper, mean, sigma, false, false, false).approximate(
            1000000L, abs_eps, -1, pedmod::cdf_methods::Korobov, 0, 8);

      expect_true(res.inform == 0L);
      expect_true(res.abserr                        < 10. * abs_eps);
      expect_true(std::abs(res.likelihood - E_prop) < 10. * abs_eps);
    }
    {
      auto res = pedmod::cdf<pedmod::likelihood>(
        func, lower, upper, mean, sigma, true, false, false).approximate(
            1000000L, abs_eps, -1, pedmod::cdf_methods::Korobov, 0, 8);

      expect_true(res.inform == 0L);
      expect_true(res.abserr                        < 10. * abs_eps);
      expect_true(std::abs(res.likelihood - E_prop) < 10. * abs_eps);
    }
    {
      auto res = pedmod::cdf<pedmod::likelihood>(
        func, lower, upper, mean, sigma, true, true, false).approximate(
            1000000L, abs_eps, -1, pedmod::cdf_methods::Korobov, 0, 8);

      expect_true(res.inform == 0L);
      expect_true(res.abserr                        < 10. * abs_eps);
      expect_true(std::abs(res.likelihood - E_prop) < 10. * abs_eps);
    }
    {
      auto res = pedmod::cdf<pedmod::likelihood>(
        func, lower, upper, mean, sigma, true, true, false).approximate(
            1000000L, abs_eps / 100, -1, pedmod::cdf_methods::Korobov, 0, 2);

      expect_true(res.inform == 0L);
      expect_true(res.abserr                        < 10. * abs_eps);
      expect_true(std::abs(res.likelihood - E_prop) < 10. * abs_eps);
    }
    {
      auto res = pedmod::cdf<pedmod::likelihood>(
        func, lower, upper, mean, sigma, true, true, false).approximate(
            100000L, abs_eps / 100, 100000L, pedmod::cdf_methods::Korobov, 0, 1);

      expect_true(res.inform == 0L);
      expect_true(res.abserr                        < 10. * abs_eps);
      expect_true(std::abs(res.likelihood - E_prop) < 10. * abs_eps);
    }
    {
      auto res = pedmod::cdf<pedmod::likelihood>(
        func, lower, upper, mean, sigma, true, true, false).approximate(
            1000000L, abs_eps, -1, pedmod::cdf_methods::Sobol, 0, 8);

      expect_true(res.inform == 0L);
      expect_true(res.abserr                        < 10. * abs_eps);
      expect_true(std::abs(res.likelihood - E_prop) < 10. * abs_eps);
    }
    {
      auto res = pedmod::cdf<pedmod::likelihood>(
        func, lower, upper, mean, sigma, true, true, false).approximate(
            1000000L, abs_eps / 100, 1000000L,
            pedmod::cdf_methods::Sobol, 0, 1);

      expect_true(res.inform == 0L);
      expect_true(std::abs(res.likelihood - E_prop) < 100. * abs_eps);
    }



    {
      auto res = pedmod::cdf<pedmod::likelihood>(
        func, lower, upper, mean, sigma, false, false, true).approximate(
            1000000L, abs_eps, -1, pedmod::cdf_methods::Korobov, 0, 8);

      expect_true(res.inform == 0L);
      expect_true(res.abserr                        < 10. * abs_eps);
      expect_true(std::abs(res.likelihood - E_prop) < 10. * abs_eps);
    }
    {
      auto res = pedmod::cdf<pedmod::likelihood>(
        func, lower, upper, mean, sigma, true, false, true).approximate(
            1000000L, abs_eps, -1, pedmod::cdf_methods::Korobov, 0, 8);

      expect_true(res.inform == 0L);
      expect_true(res.abserr                        < 10. * abs_eps);
      expect_true(std::abs(res.likelihood - E_prop) < 10. * abs_eps);
    }
    {
      auto res = pedmod::cdf<pedmod::likelihood>(
        func, lower, upper, mean, sigma, true, true, true).approximate(
            1000000L, abs_eps, -1, pedmod::cdf_methods::Korobov, 0, 8);

      expect_true(res.inform == 0L);
      expect_true(res.abserr                        < 10. * abs_eps);
      expect_true(std::abs(res.likelihood - E_prop) < 10. * abs_eps);
    }
    {
      auto res = pedmod::cdf<pedmod::likelihood>(
        func, lower, upper, mean, sigma, true, true, true).approximate(
            1000000L, abs_eps / 100, -1, pedmod::cdf_methods::Korobov, 0, 2);

      expect_true(res.inform == 0L);
      expect_true(res.abserr                        < 10. * abs_eps);
      expect_true(std::abs(res.likelihood - E_prop) < 10. * abs_eps);
    }
    {
      auto res = pedmod::cdf<pedmod::likelihood>(
        func, lower, upper, mean, sigma, true, true, true).approximate(
            100000L, abs_eps / 100, 100000L, pedmod::cdf_methods::Korobov, 0, 1);

      expect_true(res.inform == 0L);
      expect_true(res.abserr                        < 10. * abs_eps);
      expect_true(std::abs(res.likelihood - E_prop) < 10. * abs_eps);
    }
    {
      auto res = pedmod::cdf<pedmod::likelihood>(
        func, lower, upper, mean, sigma, true, true, true).approximate(
            1000000L, abs_eps, -1, pedmod::cdf_methods::Sobol, 0, 8);

      expect_true(res.inform == 0L);
      expect_true(res.abserr                        < 10. * abs_eps);
      expect_true(std::abs(res.likelihood - E_prop) < 10. * abs_eps);
    }
  }

  test_that("cdf<likelihood> gives similar output to R (1D)") {
/*
 lbs <- c(-1, -Inf, -.5)
 ubs <- c(Inf, 1, 2)
 mu <- .5
 va <- .8

 dput(mapply(function(l, u){
 f <- function(mu, va)
 pnorm(u, mean = mu, sd = sqrt(va)) - pnorm(l, mean = mu, sd = sqrt(va))
 f(mu, va)
 }, l = lbs, u = ubs))
 */
    arma::vec lbs, ubs, expect;
    constexpr double const Inf = std::numeric_limits<double>::infinity();
    lbs << -1  << -Inf << -.5;
    ubs << Inf << 1    << 2;
    expect << 0.953233743655453 << 0.711924938984711 << 0.821457505013967;
    double const mu(.5);
    double const va(.8);

    double const eps = std::pow(std::numeric_limits<double>::epsilon(), .5);
    pedmod::cdf<pedmod::likelihood>::alloc_mem(1, 1);
    pedmod::likelihood::alloc_mem(1, 1, 8);
    pedmod::likelihood func;
    for(size_t i = 0; i < 3; ++i){
      arma::vec l(1), u(1), m(1);
      arma::mat s(1, 1);
      l[0] = lbs[i];
      u[0] = ubs[i];
      m[0] = mu;
      s[0] = va;

      auto res = pedmod::cdf<pedmod::likelihood>(
        func, l, u, m, s, false, false, false).approximate(
            1000000L, 1e-8, -1, pedmod::cdf_methods::Korobov, 0, 8);

      expect_true(res.inform == 0L);
      expect_true(res.abserr                           <= 0);
      expect_true(std::abs(res.likelihood - expect[i]) <  eps);
    }

    for(size_t i = 0; i < 3; ++i){
      arma::vec l(1), u(1), m(1);
      arma::mat s(1, 1);
      l[0] = lbs[i];
      u[0] = ubs[i];
      m[0] = mu;
      s[0] = va;

      auto res = pedmod::cdf<pedmod::likelihood>(
        func, l, u, m, s, true, false, false).approximate(
            1000000L, 1e-8, -1, pedmod::cdf_methods::Korobov, 0, 8);

      expect_true(res.inform == 0L);
      expect_true(res.abserr                           <= 0);
      expect_true(std::abs(res.likelihood - expect[i]) <  eps);
    }
  }

  test_that("cdf<pedigree_l_factor> gives similar output to R (1D)") {
/*
 sigs <- list(matrix(.5, 1, 1), matrix(1.5, 1, 1))

 lws <- c(-Inf, 1  , 1)
 ubs <- c(   2, Inf, 2)
 mu <- .5
 xs <- c(.33, .5)

 f <- function(x, lw, ub, use_log){
 va <- sigs[[1]] * x[2] + sigs[[2]] * x[3] + 1
 vub <- if(is.finite(ub))
 pnorm(ub, x[1], sqrt(va), log.p = use_log) else 1
 vlb <- if(is.finite(lw))
 pnorm(lw, x[1], sqrt(va), log.p = use_log) else 0
 vub - vlb
 }

 library(numDeriv)
 dput(mapply(function(lw, ub){
 arg <- c(mu, xs)
 o <- f(arg, lw, ub, use_log = FALSE)
 do <- numDeriv::grad(f, arg, lw = lw, ub = ub, use_log = TRUE)
 c(o, do)
 }, lw = lws, ub = ubs))
 */
    arma::vec lbs, ubs;
    constexpr double const Inf = std::numeric_limits<double>::infinity();
    lbs << -Inf << 1   << 1;
    ubs << 2    << Inf << 2;
    arma::mat expect;

    expect << 0.860805198822333 << -0.186117393338058 << -0.0364459647521367 << -0.109337894247473 << 0.358932107909215 << 0.421282527277076 << 0.0274988594855953 << 0.0824965784320902 << 0.219737306731548 << 0.235165133942785 << -0.00894710527045604 << -0.0268413158128328;
    expect.reshape(4, 3);

    std::vector<arma::mat> scales;
    arma::mat s1(1, 1), s2(1, 1);
    s1.at(0, 0) =  .5;
    s2.at(0, 0) = 1.5;
    scales.emplace_back(s1);
    scales.emplace_back(s2);
    arma::mat X(1, 1, arma::fill::ones);

    pedmod::pedigree_l_factor func(scales, 1L, X, 8);
    arma::vec par;
    par << .5 << .33 << .5;

    arma::mat sig(1, 1);
    double const eps = std::pow(std::numeric_limits<double>::epsilon(), .5);
    pedmod::cdf<pedmod::pedigree_l_factor>::alloc_mem(1, 1);
    for(int i = 0; i < 3; ++i){
      func.setup(sig, par.begin() + 1, 1.);
      arma::vec lower(1), upper(1), mu(1);
      lower[0] = lbs[i];
      upper[0] = ubs[i];
      mu   [0] = par[0];

      auto const res = pedmod::cdf<pedmod::pedigree_l_factor>(
        func, lower, upper, mu, sig, true, false, false).approximate(
            1000000L, 1e-8, -1, pedmod::cdf_methods::Korobov, 0, 8);

      arma::vec const ex_res = expect.col(i);
      expect_true(res.inform == 0L);
      expect_true(res.abserr                           <= 0);

      expect_true(std::abs(res.likelihood - ex_res[0]) <  eps);
      for(int j = 0; j < 3; ++j)
        expect_true(std::abs(res.derivs[j] - ex_res[j + 1]) <  eps);
    }
  }

  test_that("cdf<pedigree_l_factor> gives similar output to R with one scale matrix") {
    /*
     sigs <- list(matrix(c(1, .25, 0, .25, 1, .1, 0, .1, 1), 3))

     lws <- c(-Inf, -1  , -1.5)
     ubs <- c(   2, Inf , 1)
     mu <- c(.5, -.25, 0)
     sc <- .5

     library(mvtnorm)
     f <- function(x){
     set.seed(1)
     mu <- x[1:3]
     Sigma <- sigs[[1L]] * x[4]
     diag(Sigma) <- diag(Sigma) + 1
     pmvnorm(lower = lws, upper = ubs, mean = mu, sigma = Sigma,
     algorithm = GenzBretz(maxpts = 1000000L, abseps = 1e-10,
     releps = 0))
     }

     dput(f(c(mu, sc)))
     library(numDeriv)
     dput(grad(f, c(mu, sc)))
     */
    std::vector<unsigned> seeds = { 1L };
    parallelrng::set_rng_seeds(seeds);

    arma::vec lbs, ubs, expect, mu;
    constexpr double const Inf = std::numeric_limits<double>::infinity();
    lbs << -Inf << -1  << -1.5;
    ubs << 2    << Inf << 1;
    expect << 0.438847591008297 << -0.0800167903883939 << 0.166356462364297 << -0.0550185761207399 << -0.186593525090286;
    mu << .5 << -.25 << 0;

    arma::mat s1;
    s1 << 1 << .25 << 0 << .25 << 1 << .1 << 0 << .1 << 1;
    s1.reshape(3, 3);
    std::vector<arma::mat> scales;
    scales.emplace_back(s1);

    arma::mat X(3, 3, arma::fill::zeros);
    X.diag().ones();
    pedmod::pedigree_l_factor func(scales, 1L, X, 8);

    arma::mat sig(3, 3);
    double const scalar = .5;
    double const eps =
      std::pow(std::numeric_limits<double>::epsilon(), .25);
    constexpr unsigned const n_deriv = 4;

    {
      // set the normalization constant
      pedmod::cdf<pedmod::likelihood>::alloc_mem(3, 1);
      pedmod::likelihood::alloc_mem(3, 1, 8);
      pedmod::likelihood lfunc;
      func.setup(sig, &scalar, 1., true);

      auto const norm_const = pedmod::cdf<pedmod::likelihood>(
        lfunc, lbs, ubs, mu, sig, false, false, false).approximate(
            1000000L, eps, -1, pedmod::cdf_methods::Korobov, 0, 8);

      func.setup(sig, &scalar, norm_const.likelihood);
    }

    pedmod::cdf<pedmod::pedigree_l_factor>::alloc_mem(3, 1);
    {
      auto const res = pedmod::cdf<pedmod::pedigree_l_factor>(
        func, lbs, ubs, mu, sig, false, false, false).approximate(
            10000000L, eps / 10000, -1, pedmod::cdf_methods::Korobov, 0, 8);

      expect_true(std::abs(res.likelihood - expect[0]) <  eps);
      expect_true(res.derivs.n_elem == n_deriv);

      for(unsigned i = 0; i < n_deriv; ++i)
        expect_true(std::abs(res.derivs[i] - expect[i + 1] / expect[0]) <  eps);
    }

    {
      auto const res = pedmod::cdf<pedmod::pedigree_l_factor>(
        func, lbs, ubs, mu, sig, true, false, false).approximate(
            10000000L, eps / 1000, -1, pedmod::cdf_methods::Korobov, 0, 8);

      expect_true(std::abs(res.likelihood - expect[0]) <  eps);
      expect_true(res.derivs.n_elem == n_deriv);
      for(unsigned i = 0; i < n_deriv; ++i)
        expect_true(std::abs(res.derivs[i] - expect[i + 1] / expect[0]) <  eps);
    }

    {
      auto const res = pedmod::cdf<pedmod::pedigree_l_factor>(
        func, lbs, ubs, mu, sig, true, true, false).approximate(
            10000000L, eps / 1000, -1, pedmod::cdf_methods::Korobov, 0, 8);

      expect_true(std::abs(res.likelihood - expect[0]) <  eps);
      expect_true(res.derivs.n_elem == n_deriv);
      for(unsigned i = 0; i < n_deriv; ++i)
        expect_true(std::abs(res.derivs[i] - expect[i + 1] / expect[0]) <  eps);
    }
    {
      auto const res = pedmod::cdf<pedmod::pedigree_l_factor>(
        func, lbs, ubs, mu, sig, true, true, false).approximate(
            10000000L, eps / 1000, -1, pedmod::cdf_methods::Korobov, 0, 2);

      expect_true(std::abs(res.likelihood - expect[0]) <  eps);
      expect_true(res.derivs.n_elem == n_deriv);
      for(unsigned i = 0; i < n_deriv; ++i)
        expect_true(std::abs(res.derivs[i] - expect[i + 1] / expect[0]) <  eps);
    }
    {
      auto const res = pedmod::cdf<pedmod::pedigree_l_factor>(
        func, lbs, ubs, mu, sig, true, true, false).approximate(
            1000000L, eps / 1000, 1000000L, pedmod::cdf_methods::Korobov, 0, 1);

      expect_true(std::abs(res.likelihood - expect[0]) <  100 * eps);
      expect_true(res.derivs.n_elem == n_deriv);
      for(unsigned i = 0; i < n_deriv; ++i)
        expect_true(std::abs(res.derivs[i] - expect[i + 1] / expect[0]) < 100 * eps);
    }
    {
      auto const res = pedmod::cdf<pedmod::pedigree_l_factor>(
        func, lbs, ubs, mu, sig, true, true, false).approximate(
            10000000L, eps / 1000, -1, pedmod::cdf_methods::Sobol, 0, 8);

      expect_true(std::abs(res.likelihood - expect[0]) <  eps);
      expect_true(res.derivs.n_elem == n_deriv);
      for(unsigned i = 0; i < n_deriv; ++i)
        expect_true(std::abs(res.derivs[i] - expect[i + 1] / expect[0]) <  eps);
    }
  }

  test_that("cdf<pedigree_l_factor> gives similar output to R with two scale matrices") {
    /*
     sigs <- list(matrix(c(1, .25, 0, .25, 1, .1, 0, .1, 1), 3),
     matrix(c(1, 0, 1, 0, 1, 0, 1, 0, 1), 3))

     lws <- c(-Inf, -1  , -1.5)
     ubs <- c(   2, Inf , 1)
     mu <- c(.5, -.25, 0)
     sc <- c(.5, .67)

     library(mvtnorm)
     f <- function(x){
     set.seed(1)
     mu <- x[1:3]
     Sigma <- sigs[[1L]] * x[4] + sigs[[2L]] * x[5]
     diag(Sigma) <- diag(Sigma) + 1
     pmvnorm(lower = lws, upper = ubs, mean = mu, sigma = Sigma,
     algorithm = GenzBretz(maxpts = 1000000L, abseps = 1e-10,
     releps = 0))
     }

     dput(f(c(mu, sc)))
     library(numDeriv)
     dput(grad(f, c(mu, sc)))
     */
    std::vector<unsigned> seeds = { 1L };
    parallelrng::set_rng_seeds(seeds);

    arma::vec lbs, ubs, expect, mu;
    constexpr double const Inf = std::numeric_limits<double>::infinity();
    lbs << -Inf << -1  << -1.5;
    ubs << 2    << Inf << 1;
    expect << 0.355656058666919 << -0.065516198469146 << 0.123908592190096 << -0.0184193612120618 << -0.116799569758003 << -0.0926881442370132;
    mu << .5 << -.25 << 0;

    arma::mat s1, s2;
    s1 << 1 << .25 << 0 << .25 << 1 << .1 << 0 << .1 << 1;
    s2 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1;
    s1.reshape(3, 3);
    s2.reshape(3, 3);
    std::vector<arma::mat> scales;
    scales.emplace_back(s1);
    scales.emplace_back(s2);

    arma::mat X(3, 3, arma::fill::zeros);
    X.diag().ones();
    pedmod::pedigree_l_factor func(scales, 1L, X, 8);

    arma::mat sig(3, 3);
    double const scs[2] = { .5, .67 };
    double const eps =
      std::pow(std::numeric_limits<double>::epsilon(), .25);
    constexpr unsigned const n_deriv = 5;

    {
      pedmod::cdf<pedmod::likelihood>::alloc_mem(3, 1);
      pedmod::likelihood::alloc_mem(3, 1, 8);
      pedmod::likelihood lfunc;
      func.setup(sig, scs, 1., true);

      auto const norm_const = pedmod::cdf<pedmod::likelihood>(
        lfunc, lbs, ubs, mu, sig, false, false, false).approximate(
            1000000L, 1e-8, -1, pedmod::cdf_methods::Korobov, 0, 8);

      func.setup(sig, scs, norm_const.likelihood);
    }

    pedmod::cdf<pedmod::pedigree_l_factor>::alloc_mem(3, 1);
    {
      auto const res = pedmod::cdf<pedmod::pedigree_l_factor>(
        func, lbs, ubs, mu, sig, false, false, false).approximate(
            10000000L, eps / 1000, -1,
            pedmod::cdf_methods::Korobov, 0, 8);

      expect_true(std::abs(res.likelihood - expect[0]) <  eps);
      expect_true(res.derivs.n_elem == n_deriv);
      for(unsigned i = 0; i < n_deriv; ++i)
        expect_true(std::abs(res.derivs[i] - expect[i + 1] / expect[0]) <  eps);
    }

    {
      auto const res = pedmod::cdf<pedmod::pedigree_l_factor>(
        func, lbs, ubs, mu, sig, true, false, false).approximate(
            10000000L, eps / 1000, -1,
            pedmod::cdf_methods::Korobov, 0, 8);

      expect_true(std::abs(res.likelihood - expect[0]) <  eps);
      expect_true(res.derivs.n_elem == n_deriv);
      for(unsigned i = 0; i < n_deriv; ++i)
        expect_true(std::abs(res.derivs[i] - expect[i + 1] / expect[0]) <  eps);
    }

    {
      auto const res = pedmod::cdf<pedmod::pedigree_l_factor>(
        func, lbs, ubs, mu, sig, true, true, false).approximate(
            10000000L, eps / 1000, -1,
            pedmod::cdf_methods::Korobov, 0, 8);

      expect_true(std::abs(res.likelihood - expect[0]) <  eps);
      expect_true(res.derivs.n_elem == n_deriv);
      for(unsigned i = 0; i < n_deriv; ++i)
        expect_true(std::abs(res.derivs[i] - expect[i + 1] / expect[0]) <  eps);
    }

    {
      auto const res = pedmod::cdf<pedmod::pedigree_l_factor>(
        func, lbs, ubs, mu, sig, true, true, false).approximate(
            10000000L, eps / 1000, -1,
            pedmod::cdf_methods::Sobol, 0, 8);

      expect_true(std::abs(res.likelihood - expect[0]) <  eps);
      expect_true(res.derivs.n_elem == n_deriv);
      for(unsigned i = 0; i < n_deriv; ++i)
        expect_true(std::abs(res.derivs[i] - expect[i + 1] / expect[0]) <  eps);
    }

    {
      auto const res = pedmod::cdf<pedmod::pedigree_l_factor>(
        func, lbs, ubs, mu, sig, true, true, false).approximate(
            10000000L, eps / 10000, -1,
            pedmod::cdf_methods::Sobol, 0, 2);

      expect_true(std::abs(res.likelihood - expect[0]) <  eps);
      expect_true(res.derivs.n_elem == n_deriv);
      for(unsigned i = 0; i < n_deriv; ++i)
        expect_true(std::abs(res.derivs[i] - expect[i + 1] / expect[0]) <  eps);
    }
  }

  test_that("cdf<generic_l_factor> gives similar output to R with a multivariate example") {
    /*
     set.seed(1)
     n <- 4
     dput(Sig <- rWishart(1, n, diag(1/n, n)) |> round(3) |> drop())
     dput(mu <- runif(n) |> round(2))

     upper_to_full <- \(x){
     dim <- (sqrt(8 * length(x) + 1) - 1) / 2
     out <- matrix(0, dim, dim)
     out[upper.tri(out, TRUE)] <- x
     out[lower.tri(out)] <- t(out)[lower.tri(out)]
     out
     }

     d_upper_to_full <- \(x){
     dim <- (sqrt(8 * length(x) + 1) - 1) / 2
     out <- matrix(0, dim, dim)
     out[upper.tri(out, TRUE)] <- x
     out[upper.tri(out)] <- out[upper.tri(out)] / 2
     out[lower.tri(out)] <- t(out)[lower.tri(out)]
     out
     }

     library(mvtnorm)
     f <- \(x){
     mu <- head(x, n)
     Sig <- tail(x, -n) |> upper_to_full()

     set.seed(1)
     pmvnorm(upper = numeric(n), mean = mu, sigma = Sig,
     algorithm = GenzBretz(maxpts = 1000000L, abseps = 0)) |>
     log()
     }

     f(c(mu, Sig[upper.tri(Sig, TRUE)])) |> exp() |> dput()
     gr <- numDeriv::grad(f, c(mu, Sig[upper.tri(Sig, TRUE)]))
     dput(c(head(gr, n), tail(gr, -n) |> d_upper_to_full()))
     */

    constexpr unsigned dim{4}, n_sequences{8};
    arma::mat Sigma{0.415, 0.41, -0.496, -0.002, 0.41, 1.791, -1.036, 1.414, -0.496, -1.036, 1.227, -0.309, -0.002, 1.414, -0.309, 1.695};
    Sigma.reshape(dim, dim);
    arma::vec const mean{0.21, 0.65, 0.13, 0.27};
    constexpr double const Inf = std::numeric_limits<double>::infinity();
    arma::vec const lb{-Inf, -Inf, -Inf, -Inf},
                    ub{0, 0, 0, 0};

    constexpr double true_likelihood{0.0069965920887278},
                     true_gr[]{-3.07835212049791, -1.56457774584047, -2.93000348516626, -0.0831394115021516, 2.95017357215864, 2.53497453851452, 3.91200497767273, 0.0400758740756828, 2.53497453851452, 0.529085992133384, 1.84325384191507, 0.304926573412647, 3.91200497767273, 1.84325384191507, 3.35128377934077, 0.23174602446223, 0.0400758740756828, 0.304926573412647, 0.23174602446223, -0.205458808913187};

    double const eps =
      std::pow(std::numeric_limits<double>::epsilon(), .25);

    std::vector<unsigned> seeds{ 1L };
    parallelrng::set_rng_seeds(seeds);

    pedmod::generic_l_factor::alloc_mem(dim, 1, n_sequences);
    pedmod::likelihood::alloc_mem(dim, 1, n_sequences);

    pedmod::cdf<pedmod::likelihood>::alloc_mem(dim, 1);
    pedmod::cdf<pedmod::generic_l_factor>::alloc_mem(dim, 1);

    pedmod::likelihood lfunc;
    auto const norm_const = pedmod::cdf<pedmod::likelihood>(
      lfunc, lb, ub, mean, Sigma, true, true, false).approximate(
          1000000L, 1e-8, -1, pedmod::cdf_methods::Korobov, 0, n_sequences);

    expect_true(
      std::abs(true_likelihood - norm_const.likelihood) <
        std::abs(true_likelihood) * eps);

    pedmod::generic_l_factor l_factor(mean, Sigma, norm_const.likelihood);
    auto const res = pedmod::cdf<pedmod::generic_l_factor>(
      l_factor, lb, ub, mean, Sigma, true, true, false).approximate(
          1000000L, -1, -1, pedmod::cdf_methods::Korobov, 0, n_sequences);

    expect_true(
      std::abs(true_likelihood - res.likelihood) <
        std::abs(true_likelihood) * eps);

    expect_true(res.derivs.size() == dim * (dim + 1));
    for(unsigned i = 0; i < dim * (dim + 1); ++i)
      expect_true(
        std::abs(res.derivs[i] - true_gr[i]) <
          std::abs(true_gr[i]) * 5e-3);
  }

  test_that("cdf<generic_l_factor> gives similar output to R with a univariate example") {
    /*
     mu <- 1.5
     Sig <- .67
     f <- \(x) pnorm(0, x[1], sqrt(x[2]), log = TRUE)
     dput(f(c(mu, Sig)) |> exp())
     dput(numDeriv::grad(f, c(mu, Sig)))
     */
    constexpr unsigned dim{4}, n_sequences{1};
    constexpr double mu{1.5},
                     Sig{.67},
                     true_likelihood{0.0334353800612323},
                     true_gr[]{-2.71919649074598, 3.04387666878185};
    constexpr double const Inf = std::numeric_limits<double>::infinity();
    arma::vec const lb{-Inf}, ub{0};

    std::vector<unsigned> seeds{ 1L };
    parallelrng::set_rng_seeds(seeds);

    pedmod::generic_l_factor::alloc_mem(dim, 1, n_sequences);
    pedmod::cdf<pedmod::generic_l_factor>::alloc_mem(dim, 1);

    arma::vec mu_vec{mu};
    arma::mat Sig_mat{Sig};
    pedmod::generic_l_factor l_factor(mu_vec, Sig_mat, 1);
    auto const res = pedmod::cdf<pedmod::generic_l_factor>(
      l_factor, lb, ub, mu_vec, Sig_mat, true, true, false).approximate(
          1000000L, -1, -1, pedmod::cdf_methods::Korobov, 0, n_sequences);

    expect_true(
      std::abs(true_likelihood - res.likelihood) <
        std::abs(true_likelihood) * 1e-8);

    expect_true(res.derivs.size() == 2);
    for(unsigned i = 0; i < 2; ++i)
      expect_true(
        std::abs(res.derivs[i] - true_gr[i]) < std::abs(true_gr[i]) * 1e-5);
  }
}
