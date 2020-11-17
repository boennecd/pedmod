#include "cdfarpx.h"
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

    double const abs_eps = std::pow(std::numeric_limits<double>::epsilon(),
                                   .33);
    double constexpr E_prop(0.0693596863013216);
    pedmod::likelihood func;
    {
      auto res = pedmod::cdf<pedmod::likelihood>(
        func, lower, upper, mean, sigma, false, false).approximate(
            1000000L, abs_eps, -1);

      expect_true(res.inform == 0L);
      expect_true(res.abserr                        < 100. * abs_eps);
      expect_true(std::abs(res.likelihood - E_prop) < 100. * abs_eps);
    }
    {
      auto res = pedmod::cdf<pedmod::likelihood>(
        func, lower, upper, mean, sigma, true, false).approximate(
            1000000L, abs_eps, -1);

      expect_true(res.inform == 0L);
      expect_true(res.abserr                        < 100. * abs_eps);
      expect_true(std::abs(res.likelihood - E_prop) < 100. * abs_eps);
    }
    {
      auto res = pedmod::cdf<pedmod::likelihood>(
        func, lower, upper, mean, sigma, true, true).approximate(
            1000000L, abs_eps, -1);

      expect_true(res.inform == 0L);
      expect_true(res.abserr                        < 100. * abs_eps);
      expect_true(std::abs(res.likelihood - E_prop) < 100. * abs_eps);
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
    pedmod::likelihood func;
    for(size_t i = 0; i < 3; ++i){
      arma::vec l(1), u(1), m(1);
      arma::mat s(1, 1);
      l[0] = lbs[i];
      u[0] = ubs[i];
      m[0] = mu;
      s[0] = va;

      auto res = pedmod::cdf<pedmod::likelihood>(
        func, l, u, m, s, false, false).approximate(1000000L, 1e-8, -1);

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
        func, l, u, m, s, true, false).approximate(1000000L, 1e-8, -1);

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

 f <- function(x, lw, ub){
 va <- sigs[[1]] * x[2] + sigs[[2]] * x[3]
 pnorm(ub, x[1], sqrt(va)) - pnorm(lw, x[1], sqrt(va))
 }

 library(numDeriv)
 dput(mapply(function(lw, ub){
 arg <- c(mu, xs)
 o <- f(arg, lw, ub)
 do <- numDeriv::grad(f, arg, lw = lw, ub = ub)
 c(o, do)
 }, lw = lws, ub = ubs))
 */
    arma::vec lbs, ubs;
    constexpr double const Inf = std::numeric_limits<double>::infinity();
    lbs << -Inf << 1   << 1;
    ubs << 2    << Inf << 2;
    arma::mat expect;

    expect << 0.941574032465999 << -0.121963784936641 << -0.0499851577637313 << -0.149955473284218 << 0.300588605141031 << 0.363805844373705 << 0.0497002519677031 << 0.149100755894646 << 0.24216263760703 << 0.241842059437064 << -0.000284905796028268 << -0.000854717389572365;
    expect.reshape(4, 3);

    std::vector<arma::mat> scales;
    arma::mat s1(1, 1), s2(1, 1);
    s1.at(0, 0) =  .5;
    s2.at(0, 0) = 1.5;
    scales.emplace_back(s1);
    scales.emplace_back(s2);

    pedmod::pedigree_l_factor func(scales);
    arma::vec par;
    par << .5 << .33 << .5;

    arma::mat sig(1, 1);
    double const eps = std::pow(std::numeric_limits<double>::epsilon(), .5);
    for(int i = 0; i < 3; ++i){
      func.setup(sig, par.begin() + 1);
      arma::vec lower(1), upper(1), mu(1);
      lower[0] = lbs[i];
      upper[0] = ubs[i];
      mu   [0] = par[0];

      auto const res = pedmod::cdf<pedmod::pedigree_l_factor>(
        func, lower, upper, mu, sig, true, false).approximate(
            1000000L, 1e-8, -1);

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
    expect << 0.757090496673361 << -0.0510313065618294 << 0.292217686124752 << -0.133769372428677 << -0.546515613536537;
    mu << .5 << -.25 << 0;

    arma::mat s1;
    s1 << 1 << .25 << 0 << .25 << 1 << .1 << 0 << .1 << 1;
    s1.reshape(3, 3);
    std::vector<arma::mat> scales;
    scales.emplace_back(s1);

    pedmod::pedigree_l_factor func(scales);

    arma::mat sig(3, 3);
    double const scalar = .5;
    func.setup(sig, &scalar);
    double const eps =
      std::pow(std::numeric_limits<double>::epsilon(), .25);
    constexpr unsigned const n_deriv = 4;

    {
      auto const res = pedmod::cdf<pedmod::pedigree_l_factor>(
        func, lbs, ubs, mu, sig, false, false).approximate(
            10000000L, eps / 10, -1);

      expect_true(std::abs(res.likelihood - expect[0]) <  eps);
      expect_true(res.derivs.n_elem == n_deriv);
      for(unsigned i = 0; i < n_deriv; ++i)
        expect_true(std::abs(res.derivs[i] - expect[i + 1]) <  eps);
    }

    {
      auto const res = pedmod::cdf<pedmod::pedigree_l_factor>(
        func, lbs, ubs, mu, sig, true, false).approximate(
            10000000L, eps / 10, -1);

      expect_true(std::abs(res.likelihood - expect[0]) <  eps);
      expect_true(res.derivs.n_elem == n_deriv);
      for(unsigned i = 0; i < n_deriv; ++i)
        expect_true(std::abs(res.derivs[i] - expect[i + 1]) <  eps);
    }

    {
      auto const res = pedmod::cdf<pedmod::pedigree_l_factor>(
        func, lbs, ubs, mu, sig, true, true).approximate(
            10000000L, eps / 10, -1);

      expect_true(std::abs(res.likelihood - expect[0]) <  eps);
      expect_true(res.derivs.n_elem == n_deriv);
      for(unsigned i = 0; i < n_deriv; ++i)
        expect_true(std::abs(res.derivs[i] - expect[i + 1]) <  eps);
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
    expect << 0.528113607673159 << -0.0635331336760654 << 0.206551124706129 << -0.0523906509478266 << -0.27157001616986 << -0.216238317754731;
    mu << .5 << -.25 << 0;

    arma::mat s1, s2;
    s1 << 1 << .25 << 0 << .25 << 1 << .1 << 0 << .1 << 1;
    s2 << 1 << 0 << 1 << 0 << 1 << 0 << 1 << 0 << 1;
    s1.reshape(3, 3);
    s2.reshape(3, 3);
    std::vector<arma::mat> scales;
    scales.emplace_back(s1);
    scales.emplace_back(s2);

    pedmod::pedigree_l_factor func(scales);

    arma::mat sig(3, 3);
    double const scs[2] = { .5, .67 };
    func.setup(sig, scs);
    double const eps =
      std::pow(std::numeric_limits<double>::epsilon(), .25);
    constexpr unsigned const n_deriv = 5;

    {
      auto const res = pedmod::cdf<pedmod::pedigree_l_factor>(
        func, lbs, ubs, mu, sig, false, false).approximate(
            10000000L, eps / 10, -1);

      expect_true(std::abs(res.likelihood - expect[0]) <  eps);
      expect_true(res.derivs.n_elem == n_deriv);
      for(unsigned i = 0; i < n_deriv; ++i)
        expect_true(std::abs(res.derivs[i] - expect[i + 1]) <  eps);
    }

    {
      auto const res = pedmod::cdf<pedmod::pedigree_l_factor>(
        func, lbs, ubs, mu, sig, true, false).approximate(
            10000000L, eps / 10, -1);

      expect_true(std::abs(res.likelihood - expect[0]) <  eps);
      expect_true(res.derivs.n_elem == n_deriv);
      for(unsigned i = 0; i < n_deriv; ++i)
        expect_true(std::abs(res.derivs[i] - expect[i + 1]) <  eps);
    }

    {
      auto const res = pedmod::cdf<pedmod::pedigree_l_factor>(
        func, lbs, ubs, mu, sig, true, true).approximate(
            10000000L, eps / 10, -1);

      expect_true(std::abs(res.likelihood - expect[0]) <  eps);
      expect_true(res.derivs.n_elem == n_deriv);
      for(unsigned i = 0; i < n_deriv; ++i)
        expect_true(std::abs(res.derivs[i] - expect[i + 1]) <  eps);
    }
  }
}
