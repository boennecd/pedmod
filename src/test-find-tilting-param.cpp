#include <testthat.h>
#include "find-tilting-param.h"
#include <limits.h>

namespace {
constexpr double Inf{std::numeric_limits<double>::infinity()};
}

context("find_tilting_param unit tests") {
  test_that("find_tilting_param gives the same as an R implementation in a 5D example") {
    /*
     .check_input <- \(n, a, b, Sig)
     stopifnot(length(a) == n, length(b) == n, all(dim(Sig) == c(n, n)),
     all(a < b))

     psi <- \(x, a, b, Sig){
     n <- length(x) %/% 2L
     mu <- head(x, n)
     x <- tail(x, n)

     .check_input(n, a, b, Sig)

     ubs <- lbs <- numeric(n)
     C <- chol(Sig)
     lbs[1] <- a[1] / C[1, 1]
     ubs[1] <- b[1] / C[1, 1]

     for(i in seq_len(n - 1L) + 1L){
     lbs[i] <- (a[i] - sum(C[seq_len(i - 1L), i] * x[seq_len(i - 1L)])) / C[i, i]
     ubs[i] <- (b[i] - sum(C[seq_len(i - 1L), i] * x[seq_len(i - 1L)])) / C[i, i]
     }
     lbs[is.infinite(a)] <- -Inf
     ubs[is.infinite(b)] <-  Inf

     -sum(mu * x) + sum(mu^2) / 2 + sum(log(pnorm(ubs - mu) - pnorm(lbs - mu)))
     }

     psi_safe <- \(x, a, b, Sig){
     n <- length(x) %/% 2L
     mu <- head(x, n)
     x <- tail(x, n)

     .check_input(n, a, b, Sig)

     ubs <- lbs <- numeric(n)
     C <- chol(Sig)
     lbs[1] <- a[1] / C[1, 1]
     ubs[1] <- b[1] / C[1, 1]

     for(i in seq_len(n - 1L) + 1L){
     lbs[i] <- (a[i] - sum(C[seq_len(i - 1L), i] * x[seq_len(i - 1L)])) / C[i, i]
     ubs[i] <- (b[i] - sum(C[seq_len(i - 1L), i] * x[seq_len(i - 1L)])) / C[i, i]
     }

     pnrm_log_lbs <- ifelse(is.infinite(lbs), -Inf, pnorm(lbs - mu, log.p = TRUE))
     pnrm_log_ubs <- ifelse(is.infinite(ubs), Inf, pnorm(ubs - mu, log.p = TRUE))

     pnrm_terms <- ifelse(
     is.infinite(ubs),
     log1p(-exp(pnrm_log_lbs)),
     pnrm_log_ubs + log1p(-exp(pnrm_log_lbs - pnrm_log_ubs)))

     -sum(mu * x) + sum(mu^2) / 2 + sum(pnrm_terms)
     }

     d_psi <- \(x, a, b, Sig){
     n <- length(x) %/% 2L
     mu <- head(x, n)
     x <- tail(x, n)

     .check_input(n, a, b, Sig)

     ubs <- lbs <- numeric(n)
     C <- chol(Sig)
     lbs[1] <- a[1] / C[1, 1]
     ubs[1] <- b[1] / C[1, 1]

     for(i in seq_len(n - 1L) + 1L){
     lbs[i] <- (a[i] - sum(C[seq_len(i - 1L), i] * x[seq_len(i - 1L)])) / C[i, i]
     ubs[i] <- (b[i] - sum(C[seq_len(i - 1L), i] * x[seq_len(i - 1L)])) / C[i, i]
     }

     pnrm_log_lbs <- ifelse(is.infinite(lbs), -Inf, pnorm(lbs - mu, log.p = TRUE))
     pnrm_log_ubs <- ifelse(is.infinite(ubs), Inf, pnorm(ubs - mu, log.p = TRUE))

     denoms_log <- ifelse(
     is.infinite(ubs),
     log1p(-exp(pnrm_log_lbs)),
     pnrm_log_ubs + log1p(-exp(pnrm_log_lbs - pnrm_log_ubs)))

     dnrms_log_lbs <- dnorm(lbs - mu, log = TRUE)
     dnrms_log_ubs <- dnorm(ubs - mu, log = TRUE)

     ratio_lbs <- exp(dnrms_log_lbs - denoms_log)
     ratio_ubs <- exp(dnrms_log_ubs - denoms_log)

     derivs <- ratio_lbs - ratio_ubs
     C <- C %*% diag(diag(C)^-1)

     c(mu - x + derivs, -mu + (C - diag(n)) %*% derivs)
     }

     dd_psi <- \(x, a, b, Sig){
     n <- length(x) %/% 2L
     mu <- head(x, n)
     x <- tail(x, n)

     .check_input(n, a, b, Sig)

     ubs <- lbs <- numeric(n)
     C <- chol(Sig)
     lbs[1] <- a[1] / C[1, 1]
     ubs[1] <- b[1] / C[1, 1]

     for(i in seq_len(n - 1L) + 1L){
     lbs[i] <- (a[i] - sum(C[seq_len(i - 1L), i] * x[seq_len(i - 1L)])) / C[i, i]
     ubs[i] <- (b[i] - sum(C[seq_len(i - 1L), i] * x[seq_len(i - 1L)])) / C[i, i]
     }

     pnrm_log_lbs <- ifelse(is.infinite(lbs), -Inf, pnorm(lbs - mu, log.p = TRUE))
     pnrm_log_ubs <- ifelse(is.infinite(ubs), Inf, pnorm(ubs - mu, log.p = TRUE))

     denoms_log <- ifelse(
     is.infinite(ubs),
     log1p(-exp(pnrm_log_lbs)),
     pnrm_log_ubs + log1p(-exp(pnrm_log_lbs - pnrm_log_ubs)))

     dnrms_log_lbs <- dnorm(lbs - mu, log = TRUE)
     dnrms_log_ubs <- dnorm(ubs - mu, log = TRUE)

     ratio_lbs <- exp(dnrms_log_lbs - denoms_log)
     ratio_ubs <- exp(dnrms_log_ubs - denoms_log)

     derivs_gr <- ratio_lbs - ratio_ubs

     derivs <-
     ifelse(is.infinite(lbs), 0, (lbs - mu) * ratio_lbs) -
     ifelse(is.infinite(ubs), 0, (ubs - mu) * ratio_ubs) -
     derivs_gr^2
     C <- C %*% diag(diag(C)^-1)
     diff_mat <- C - diag(n)

     out <- matrix(0., 2L * n, 2L * n)
     m1 <- diff_mat %*% diag(derivs)
     out[  1:n ,   1:n ] <- diag(n) + diag(derivs)
     out[-(1:n),   1:n ] <- m1 - diag(n)
     out[-(1:n), -(1:n)] <- tcrossprod(m1, diff_mat)
     out[upper.tri(out)] <- t(out)[upper.tri(out)]
     out
     }

# examples
     set.seed(111)
     n <- 5
     Sig <- matrix(1, n, n)
     while (cov2cor(Sig)[lower.tri(Sig)] |> abs() |> max() > .999)
     Sig <- rWishart(1, n, diag(n)) |> drop()
     a <- runif(n, -2, 0)
     b <- runif(n, 0, 2)
     type <- sample.int(3, n, replace = TRUE)
     a[type == 1] <- -Inf
     b[type == 2] <- Inf

     start <- local({
     C <- chol(Sig)
     start_org_scale <- (cbind(a, b) |> rowMeans()) / diag(C)
     start_org_scale[type == 1] <- b[type == 1] / diag(C)[type == 1] - 1
     start_org_scale[type == 2] <- a[type == 2] / diag(C)[type == 2] + 1
     start_org_scale <- start_org_scale * diag(C)
     solve(t(C), start_org_scale)
     })

     par <- local({
     mu <- numeric(n)
     C <- chol(Sig)
     for(i in seq_len(n - 1L) + 1L)
     mu[i] <- -C[seq_len(i - 1L), i] %*% start[seq_len(i - 1L)] / C[i, i]
     c(mu, start)
     })

# do we start of in an interior point?
     ptr <- crossprod(chol(Sig), tail(par, n))
     all(ptr > a)
     all(ptr < b)

     psi(par, a, b, Sig)
     all.equal(psi(par, a, b, Sig), psi_safe(par, a, b, Sig))
     psi <- psi_safe

     stopifnot(all.equal(d_psi(par, a, b, Sig),
     numDeriv::grad(psi, par, a = a, b = b, Sig = Sig)))

     stopifnot(all.equal(dd_psi(par, a, b, Sig),
     numDeriv::jacobian(d_psi, par, a = a, b = b, Sig = Sig)))

# finds a root
     root_finder <- \(x, a, b, Sig, abstol = 1e-2){
     f <- \(x) d_psi(x, a, b, Sig)^2 |> sum()
     d_f <- \(x){
     f_vals <- d_psi(x, a, b, Sig)
     grs <- dd_psi(x, a, b, Sig)
     2 * rowSums(grs %*% diag(f_vals))
     }

# sanity check as this is still experimental
     num_gr <- try(numDeriv::grad(f, x), silent = TRUE)
     if(!inherits(num_gr, "try-error")){
     is_equal <- all.equal(d_f(x), num_gr, tolerance = 1e-5)
     if(!isTRUE(is_equal))
     warning(paste0(capture.output(is_equal), collapse = "\n"))
     }

# find the root
     optim(x, fn = f, gr = d_f, method = "BFGS",
 control = list(reltol = 1e-8, abstol = abstol))
 }

 res <- root_finder(par, a, b, Sig, 0)
 rbind(Estimate = res$par, Start = par)
 d_psi(res$par, a, b, Sig) |> abs() |> sum() # ~ zero
 res$counts

# do we have an interior solution
 ptr <- crossprod(chol(Sig), tail(res$par, n))
 all(ptr > a)
 all(ptr < b)

 head(res$par, n) |> dput()
 dput(a)
 dput(b)
 C <- chol(Sig)
 dput((C %*% diag(diag(C)^-1))[upper.tri(C, TRUE)])
     */
    constexpr size_t dim{5};
    constexpr double lower_limits[]{-0.815837000441616, -0.288358708962961, -Inf, -Inf, -Inf},
                     upper_limits[]{Inf, Inf, 1.26198812572033, 0.42491765901508, 0.670519533487418},
                         cholesky[]{1, -0.137599036720516, 1, -0.895959210162758, 0.155032252997349, 1, -0.0641444623930221, -0.150172397352036, 0.681663684155795, 1, 0.729181873459802, -1.43240521262587, -0.0784937458067778, -0.328361932329834, 1},
                             tilt[]{-0.070383996090136, 0.273675444869523, -0.252692930348561, 0.0529469150197525, -1.66721318828181e-11};

    auto res =  find_tilting_param
      (dim, lower_limits, upper_limits, cholesky, 1e-8);

    expect_true(res.success);
    expect_true(res.tilting_param.size() == dim);
    constexpr double eps{1e-5};
    for(size_t i = 0; i < dim; ++i)
      expect_true
        (std::abs(res.tilting_param[i] - tilt[i]) <
          (eps + std::abs(tilt[i])) * eps);
  }
}
