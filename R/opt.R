#' Optimize the Log Marginal Likelihood
#'
#' Optimizes \code{\link{eval_pedigree_ll}} and \code{\link{eval_pedigree_grad}}
#' using a passed optimization function.
#'
#' \code{pedmod_start} yield starting values which can be used for
#' \code{pedmod_opt}. The method is based on a heuristic where we make the
#' assumption that the fixed effects are unrelated to the random effects.
#'
#' @inheritParams eval_pedigree_ll
#' @param par starting values passed to \code{opt_func}.
#' @param opt_func function to perform minimization with arguments like
#' \code{\link{optim}}. BFGS is used with \code{\link{optim}} if this argument
#' is \code{NULL}.
#' @param seed seed to pass to \code{\link{set.seed}} before each gradient and
#' function evaluation. Use \code{NULL} if the seed should not be fixed.
#' @param fix integer vector with indices of \code{par} to fix. This is useful
#' for computing profile likelihoods. \code{NULL} yields all parameters.
#' @param ... Arguments passed to \code{opt_func}.
#' @param abs_eps absolute convergence threshold for
#' \code{\link{eval_pedigree_ll}} and \code{\link{eval_pedigree_grad}}.
#' @param rel_eps rel_eps convergence threshold for
#' \code{\link{eval_pedigree_ll}} and \code{\link{eval_pedigree_grad}}.
#'
#' @seealso
#' \code{\link{pedmod_sqn}} and \code{\link{pedmod_start}}.
#'
#' @return
#' \code{pedmod_opt}: The output from the \code{opt_func} argument. Thus, if
#' \code{fix} is supplied then this is optimal values of only \code{par[-fix]}
#' with
#' \code{par[fix]} being fixed to the inputs. Thus, the length is only the
#' number of non-fixed parameters.
#'
#' @importFrom stats optim
#' @export
pedmod_opt <- function(ptr, par, maxvls, abs_eps, rel_eps,
                       opt_func = NULL, seed = 1L, indices = NULL, minvls = -1L,
                       do_reorder = TRUE, use_aprx = FALSE, n_threads = 1L,
                       cluster_weights = NULL, fix = NULL, standardized = FALSE,
                       ...){
  # handle defaults
  if(is.null(opt_func)){
    opt_func <- optim
    formals(opt_func)$method <- "BFGS"
  }
  if(is.null(fix))
    fix <- integer()
  else
    fix <- unique(fix)

  # checks
  stopifnot(is.function(opt_func),
            is.integer(fix), all(fix %in% seq_along(par)),
            is.numeric(par))

  # assign objective and gradient function
  get_par <- function(x){
    if(length(fix) < 1L)
      return(x)
    out <- par
    out[-fix] <- x
    out
  }
  fn <- function(x){
    if(!is.null(seed))
      set.seed(seed)
    out <- try(-eval_pedigree_ll(
        ptr = ptr, par = get_par(x), maxvls = maxvls, rel_eps = rel_eps,
        indices = indices, minvls = minvls, abs_eps = abs_eps,
        do_reorder = do_reorder, use_aprx = use_aprx, n_threads = n_threads,
        cluster_weights = cluster_weights, standardized = standardized),
      silent = TRUE)
    if(inherits(out, "try-error"))
      return(NA_real_)
    out
  }
  gr <- function(x){
    if(!is.null(seed))
      set.seed(seed)
    out <-
      -eval_pedigree_grad(ptr = ptr, par = get_par(x), maxvls = maxvls,
                          rel_eps = rel_eps, indices = indices, minvls = minvls,
                          abs_eps = abs_eps, do_reorder = do_reorder,
                          use_aprx = use_aprx, n_threads = n_threads,
                          cluster_weights = cluster_weights,
                          standardized = standardized)
    if(length(fix) > 0)
      out <- out[-fix]
    out
  }

  # optimize and return
  opt_func(par = if(length(fix) > 0) par[-fix] else par, fn = fn, gr = gr, ...)
}

#' @rdname pedmod_opt
#' @param data the \code{\link{list}} that was passed to
#' \code{\link{get_pedigree_ll_terms}}.
#' @param scale_max the maximum value for the scale parameters. Sometimes, the
#' optimization method tends to find large scale parameters and get stuck.
#' Setting a maximum solves this.
#'
#' @return
#' \code{pedmod_start}: A \code{list} with:
#' \itemize{
#'  \item par: the starting value.
#'  \item beta_no_rng: the fixed effects MLEs without random effects.
#'  \item logLik_no_rng: the log maximum likelihood without random effects.
#'  \item logLik_est: the likelihood at par.
#' }
#' @importFrom stats binomial glm.fit
#' @importFrom utils tail
#' @export
pedmod_start <- function(ptr, data, maxvls = 1000L, abs_eps = 0, rel_eps = 1e-2,
                         seed = 1L, indices = NULL, scale_max = 9,
                         minvls = 100L, do_reorder = TRUE, use_aprx = TRUE,
                         n_threads = 1L, cluster_weights = NULL,
                         standardized = FALSE){
  # checks
  stopifnot(is.numeric(scale_max), length(scale_max) == 1L,
            scale_max > 0,
            is.numeric(seed), length(seed) %in% 0:1)

  #####
  # get the starting values for the fixed effects
  y <- unlist(lapply(data, `[[`, "y"))
  X <- do.call(rbind, lapply(data, `[[`, "X"))
  if(!is.null(cluster_weights))
    w <- unlist(Map(
      rep, cluster_weights,
      times = sapply(data, function(x) length(x$y))))
  else
    w <- rep(1, length(y))

  # checks
  n <- length(y)
  stopifnot(
    length(y) > 0, is.numeric(y),
    NROW(X) == n,
    is.null(cluster_weights) || (
      length(cluster_weights) == length(data) &&
        is.numeric(cluster_weights)))

  # fit the model without random effects
  start_fit <-  glm.fit(X, y, weights = w, family = binomial("probit"))
  beta <- start_fit$coefficients
  logLik_no_rng <- -sum(start_fit$deviance) / 2

  if(standardized){
    fn <- function(sc){
      if(!is.null(seed))
        set.seed(seed)

      sc_direct <- standardized_to_direct(sc, length(sc))
      if(any(sc_direct > scale_max))
        return(NA_real_)

      out <- try(-eval_pedigree_ll(
        ptr = ptr, par = c(beta, sc), maxvls = maxvls, abs_eps = abs_eps,
        rel_eps = rel_eps, minvls = minvls, use_aprx = use_aprx,
        n_threads = n_threads, cluster_weights = cluster_weights,
        indices = indices, do_reorder = do_reorder,
        standardized = standardized),
        silent = TRUE)
      if(inherits(out, "try-error"))
        return(NA_real_)

      out
    }
    gr <- function(sc){
      if(!is.null(seed))
        set.seed(seed)

      out <- -eval_pedigree_grad(
        ptr = ptr, par = c(beta, sc), maxvls = maxvls, abs_eps = abs_eps,
        rel_eps = rel_eps, minvls = minvls, use_aprx = use_aprx,
        n_threads = n_threads, cluster_weights = cluster_weights,
        indices = indices, do_reorder = do_reorder,
        standardized = standardized)

      tail(out, -length(beta))
    }

    # optimize
    n_scales <- get_n_scales(ptr)
    stopifnot(n_scales > 0)

    for(sc_sqrt in seq(.1, sqrt(scale_max), length.out = 5)){
      sc <- rep(log(sc_sqrt) * 2, n_scales)
      opt <- try(optim(sc, fn, gr, method = "BFGS",
                       control = list(maxit = 1000L)),
                 silent = TRUE)
      if(!inherits(opt, "try-error") && opt$convergence < 2)
        # found a good solution
        break
    }

    if(inherits(opt, "try-error") || opt$convergence > 1){
      # fall to find a good solution
      sc <- rep(log(.01), n_scales)
      logLik_est <- -fn(sc)

    } else {
      sc <- opt$par
      logLik_est <- -opt$value

    }

    return(list(par = c(beta, sc), beta_no_rng = beta,
                logLik_no_rng = logLik_no_rng,
                logLik_est = logLik_est, opt = opt))
  }

  #####
  # optimize the model where beta is proportional to estimates without the
  # random effects
  fn <- function(sc){
    if(!is.null(seed))
      set.seed(seed)
    beta_scaled <- beta * sqrt(1 + sum(exp(sc)))

    out <- try(-eval_pedigree_ll(
      ptr = ptr, par = c(beta_scaled, sc), maxvls = maxvls, abs_eps = abs_eps,
      rel_eps = rel_eps, minvls = minvls, use_aprx = use_aprx,
      n_threads = n_threads, cluster_weights = cluster_weights,
      indices = indices, do_reorder = do_reorder),
      silent = TRUE)
    if(inherits(out, "try-error"))
      return(NA_real_)

    out
  }
  gr <- function(sc){
    if(!is.null(seed))
      set.seed(seed)
    fac <- sqrt(1 + sum(exp(sc)))
    beta_scaled <- beta * fac

    out <- -eval_pedigree_grad(
      ptr = ptr, par = c(beta_scaled, sc), maxvls = maxvls, abs_eps = abs_eps,
      rel_eps = rel_eps, minvls = minvls, use_aprx = use_aprx,
      n_threads = n_threads, cluster_weights = cluster_weights,
      indices = indices, do_reorder = do_reorder)

    sum_d_beta <- sum(beta * out[seq_along(beta)])
    sum_d_beta * exp(sc) / (2 * fac) + tail(out, -length(beta))
  }

  # optimize
  n_scales <- get_n_scales(ptr)
  stopifnot(n_scales > 0)

  for(sc_sqrt in seq(.1, sqrt(scale_max), length.out = 5)){
    sc <- rep(2 * log(sc_sqrt), n_scales)
    opt <- try(optim(sc, fn, gr, upper = rep(log(scale_max), n_scales),
                     method = "L-BFGS-B",
                     control = list(lmm = 10L, maxit = 1000L)),
               silent = TRUE)
    if(!inherits(opt, "try-error") && opt$convergence < 2)
      # found a good solution
      break
  }

  if(inherits(opt, "try-error") || opt$convergence > 1){
    # fall to find a good solution
    sc <- rep(log(.01), n_scales)
    beta_scaled <- beta * sqrt(1 + sum(exp(sc)))
    logLik_est <- -fn(sc)

  } else {
    sc <- opt$par
    beta_scaled <- beta * sqrt(1 + sum(exp(sc)))
    logLik_est <- -opt$value

  }

  # return
  list(par = c(beta_scaled, sc), beta_no_rng = beta,
       logLik_no_rng = logLik_no_rng,
       logLik_est = logLik_est, opt = opt)
}

#' Optimize the Log Marginal Likelihood Using a Stochastic Quasi-Newton Method
#'
#' Optimizes \code{\link{eval_pedigree_ll}} and \code{\link{eval_pedigree_grad}}
#' using a stochastic quasi-Newton method.
#'
#' @inheritParams pedmod_opt
#' @param par starting values.
#' @param step_factor factor used for the step size. The step size is
#' \code{step_factor} divided by the iteration number.
#' @param n_it number of stochastic gradient steps to make.
#' @param n_grad_steps number of stochastic gradient steps to make between each
#' Hessian approximation update.
#' @param n_grad number of log marginal likelihood terms to include in the
#' stochastic gradient step.
#' @param n_hess number of log marginal likelihood terms to include in the
#' gradients used for the Hessian approximation update. This is set to the
#' entire sample (or \code{indices}) if this is greater than or equal to half
#' the number of log marginal likelihood terms.
#' @param minvls_hess \code{minvls} argument to use when updating the Hessian
#' approximation.
#' @param maxvls_hess \code{maxvls} argument to use when updating the Hessian
#' approximation.
#' @param abs_eps_hess \code{abs_eps} argument to use when updating the Hessian
#' approximation.
#' @param rel_eps_hess \code{rel_eps} argument to use when updating the Hessian
#' approximation.
#' @param verbose logical for whether to print output during the estimation.
#'
#' @details
#' The function uses a stochastic quasi-Newton method like suggested by
#' Byrd et al. (2016) with a few differences: Differences in gradients are
#' used rather than Hessian-vector products, BFGS rather than L-BFGS is used
#' because the problem is typically low dimensional, and damped BFGS updates
#' are used (see e.g. chapter 18 of Nocedal and Wright, 2006).
#'
#' Separate arguments for the gradient approximation in the Hessian update are
#' provided as one may want a more precise approximation for these gradients.
#' \code{step_factor} likely depends on the other parameters and the data set
#' and should be altered.
#'
#' @references
#'
#' Byrd, R. H., Hansen, S. L., Nocedal, J., & Singer, Y. (2016).
#' \emph{A stochastic quasi-Newton method for large-scale optimization}.
#' SIAM Journal on Optimization, 26(2), 1008-1031.
#'
#' Nocedal, J., & Wright, S. (2006). \emph{Numerical optimization}.
#' Springer Science & Business Media.
#'
#' @seealso
#' \code{\link{pedmod_opt}} and \code{\link{pedmod_start}}.
#'
#' @export
pedmod_sqn <- function(ptr, par, maxvls, abs_eps, rel_eps, step_factor,
                       n_it, n_grad_steps, indices = NULL, minvls = -1L,
                       n_grad = 50L, n_hess = 500L, do_reorder = TRUE,
                       use_aprx = FALSE, n_threads = 1L, cluster_weights = NULL,
                       fix = NULL, standardized = FALSE, minvls_hess = minvls,
                       maxvls_hess = maxvls, abs_eps_hess = abs_eps,
                       rel_eps_hess = rel_eps, verbose = FALSE){
  #####
  # setup before the estimation
  n_pars <- length(par) - length(fix)
  omegas <- matrix(NA_real_, n_pars, floor(n_it / n_grad_steps) + 1L)
  H <- diag(n_pars)
  any_fixed <- length(fix) > 0
  w_old <- if(any_fixed) par[-fix] else par
  omegas[, 1L] <- w_old

  n_terms <- get_n_terms(ptr)
  if(is.null(indices))
    indices <- 0:(n_terms - 1L)
  n_terms <- min(n_terms, length(indices))
  n_grad <- min(n_grad, n_terms)
  n_hess <- min(n_hess, n_terms)

  # we alter n_hess if it is greater than or equal to the number of observations
  if(n_hess >= n_terms / 2)
    n_hess <- n_terms
  hess_use_all <- n_hess == n_terms

  # assign the function and the gradient functions
  get_par <- function(x){
    if(!any_fixed)
      return(x)
    out <- par
    out[-fix] <- x
    out
  }
  fn <- function(x, abs_eps, rel_eps, minvls, maxvls, indices){
    out <- try(-eval_pedigree_ll(
      ptr = ptr, par = get_par(x), maxvls = maxvls, rel_eps = rel_eps,
      indices = indices, minvls = minvls, abs_eps = abs_eps,
      do_reorder = do_reorder, use_aprx = use_aprx, n_threads = n_threads,
      cluster_weights = cluster_weights, standardized = standardized),
      silent = TRUE)
    if(inherits(out, "try-error"))
      return(NA_real_)
    if(!is.null(out))
      out / length(indices) else out / n_terms
  }
  gr <- function(x, abs_eps, rel_eps, minvls, maxvls, indices){
    out <-
      -eval_pedigree_grad(ptr = ptr, par = get_par(x), maxvls = maxvls,
                          rel_eps = rel_eps, indices = indices, minvls = minvls,
                          abs_eps = abs_eps, do_reorder = do_reorder,
                          use_aprx = use_aprx, n_threads = n_threads,
                          cluster_weights = cluster_weights,
                          standardized = standardized)
    if(any_fixed)
      out <- out[-fix]
    if(!is.null(indices))
      out / length(indices) else out / n_terms
  }

  if(hess_use_all)
    # we will need this later
    g_new <- gr(omegas[, 1], abs_eps = abs_eps_hess, rel_eps = rel_eps_hess,
                minvls = minvls_hess, maxvls = maxvls_hess,
                indices = indices)

  #####
  # perform the optimization
  k <- 0L
  stopifnot(n_it > 0, n_grad_steps > 0)
  w_vals <- matrix(NA_real_, n_pars, n_grad_steps)
  w_new <- w_old
  t <- 1L
  while(k < n_it){
    w_vals[] <- NA_real_
    for(i in 1:n_grad_steps){
      if((k <- k + 1L) > n_it)
        break

      # perform the gradient step
      w_old <- w_new
      S <- sample(indices, n_grad)
      gr_val <- gr(w_old, abs_eps = abs_eps, rel_eps = rel_eps,
                   minvls = minvls, maxvls = maxvls, indices = S)
      w_new <- w_old - step_factor / k * drop(H %*% gr_val)
      w_vals[, i] <- w_new
    }

    if(i < n_grad_steps)
      # no need for an update of the Hessian approximation
      break

    t <- t + 1L
    omegas[, t] <- rowMeans(w_vals)
    s <- omegas[, t] - omegas[, t - 1L]
    S <- if(hess_use_all) indices else sample(indices, n_hess)

    . <- function(x)
      gr(x, abs_eps = abs_eps_hess, rel_eps = rel_eps_hess,
         minvls = minvls_hess, maxvls = maxvls_hess, indices = S)
    if(hess_use_all){
      g_old <- g_new
      g_new <- .(omegas[, t     ])

    } else {
      g_old <- .(omegas[, t - 1L])
      g_new <- .(omegas[, t     ])
    }

    if(verbose){
      cat(sprintf(
        "\nHessian update %5d\nLog marignal likelihood at previous and new average parameter vector on the same sample is:\n  %14.3f\n  %14.3f\n",
        t, attr(g_old, "logLik"), attr(g_new, "logLik")))
      cat("New average parameter vector is\n")
      print(omegas[, t])
    }

    y <- g_new - g_old
    s_y <- sum(s  * y)
    if(t <= 2L){
      # the first iteration
      rho <- s_y / sum(y * y)
      if(is.finite(rho) && rho > 0)
        H <- diag(rho, n_pars)
    }

    # damped BFGS update
    B_s <- solve(H, s)
    s_B_s <- drop(s %*% B_s)
    theta <- if(s_y >= .2 * s_B_s)
      1 else .8 * s_B_s / (s_B_s - s_y)
    r <- theta * y + (1 - theta) * B_s

    # TODO: can be done smarter
    r_s <- sum(r * s)
    D <- -outer(r, s) / r_s
    diag(D) <- diag(D) + 1

    H <- crossprod(D, H %*% D) + outer(s, s) / r_s
  }

  list(par = w_new, omegas = omegas, H = H)
}
