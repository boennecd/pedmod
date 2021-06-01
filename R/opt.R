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
                       method = 0L, ...){
  # handle defaults
  if(is.null(opt_func)){
    opt_func <- optim
    formals(opt_func)$method <- "BFGS"
    formals(opt_func)$control <- list(fnscale = get_n_terms(ptr))
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
        cluster_weights = cluster_weights, standardized = standardized,
        method = method),
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
                          standardized = standardized,
                          method = method)
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
#' @param sc_start starting value for the scale parameters. Use \code{NULL} if
#' you have no value to start with.
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
                         standardized = FALSE, method = 0L,
                         sc_start = NULL){
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

  # function to optimize the parameters
  do_opt <- function(sc, fn, gr){
    n_scales <- get_n_scales(ptr)
    optim(sc, fn, gr, upper = rep(log(scale_max), n_scales),
          method = "L-BFGS-B", control = list(
            lmm = 10L, maxit = 1000L, fnscale = get_n_terms(ptr)))
  }
  try_opt <- function(fn, gr){
    n_scales <- get_n_scales(ptr)
    stopifnot(n_scales > 0,
              is.null(sc_start) ||
                (length(sc_start) == n_scales && is.numeric(sc_start) &&
                   all(is.finite(sc_start) & sc_start > 0 &
                         sc_start <= scale_max)))

    check_ok <- function(opt)
      !inherits(opt, "try-error") && opt$convergence < 2

    has_starting_value <- !is.null(sc_start)
    if(has_starting_value){
      sc <- log(sc_start)
      opt <- try(do_opt(sc, fn = fn, gr = gr), silent = TRUE)
      is_ok <- check_ok(opt)
    }

    if(!has_starting_value || !is_ok)
      for(sc_sqrt in seq(.1, sqrt(scale_max), length.out = 5)){
        sc <- rep(log(sc_sqrt) * 2, n_scales)
        opt <- try(do_opt(sc, fn = fn, gr = gr), silent = TRUE)
        if(is_ok <- check_ok(opt))
          # found a good solution
          break
      }

    list(opt = opt, sc = sc, is_ok = is_ok)
  }


  if(standardized){
    fn <- function(sc){
      if(!is.null(seed))
        set.seed(seed)

      out <- try(-eval_pedigree_ll(
        ptr = ptr, par = c(beta, sc), maxvls = maxvls, abs_eps = abs_eps,
        rel_eps = rel_eps, minvls = minvls, use_aprx = use_aprx,
        n_threads = n_threads, cluster_weights = cluster_weights,
        indices = indices, do_reorder = do_reorder,
        standardized = standardized, method = method),
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
        standardized = standardized, method = method)

      tail(out, -length(beta))
    }

    # optimize
    res <- try_opt(fn = fn, gr = gr)

    if(!res$is_ok){
      # fall to find a good solution
      n_scales <- get_n_scales(ptr)
      sc <- rep(log(.01), n_scales)
      logLik_est <- -fn(sc)

    } else {
      sc <- res$opt$par
      logLik_est <- -res$opt$value
    }

    return(list(par = c(beta, sc), beta_no_rng = beta,
                logLik_no_rng = logLik_no_rng,
                logLik_est = logLik_est, opt = res$opt))
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
      indices = indices, do_reorder = do_reorder, method = method),
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
      indices = indices, do_reorder = do_reorder, method = method)

    sum_d_beta <- sum(beta * out[seq_along(beta)])
    sum_d_beta * exp(sc) / (2 * fac) + tail(out, -length(beta))
  }

  # optimize
  res <- try_opt(fn = fn, gr = gr)

  if(!res$is_ok){
    # fall to find a good solution
    n_scales <- get_n_scales(ptr)
    sc <- rep(log(.01), n_scales)
    beta_scaled <- beta * sqrt(1 + sum(exp(sc)))
    logLik_est <- -fn(sc)

  } else {
    sc <- res$opt$par
    beta_scaled <- beta * sqrt(1 + sum(exp(sc)))
    logLik_est <- -res$opt$value

  }

  # return
  list(par = c(beta_scaled, sc), beta_no_rng = beta,
       logLik_no_rng = logLik_no_rng,
       logLik_est = logLik_est, opt = res$opt)
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
#' @param check_every integer for the number of gradient steps between checking
#' that the likelihood did increase. If not, the iterations are reset and the
#' step-size is halved.
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
                       rel_eps_hess = rel_eps, verbose = FALSE,
                       method = 0L, check_every = 2L * n_grad_steps){
  #####
  # setup before the estimation
  n_pars <- length(par) - length(fix)
  omegas <- matrix(NA_real_, n_pars, ceiling(n_it / n_grad_steps) + 2L)
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
      cluster_weights = cluster_weights, standardized = standardized,
      method = method),
      silent = TRUE)
    if(inherits(out, "try-error"))
      return(NA_real_)

    denom <- if(!is.null(indices))
      length(indices) else n_terms

    structure(c(out) / denom, std = attr(out, "std") / denom)
  }
  gr <- function(x, abs_eps, rel_eps, minvls, maxvls, indices){
    out <-
      -eval_pedigree_grad(ptr = ptr, par = get_par(x), maxvls = maxvls,
                          rel_eps = rel_eps, indices = indices, minvls = minvls,
                          abs_eps = abs_eps, do_reorder = do_reorder,
                          use_aprx = use_aprx, n_threads = n_threads,
                          cluster_weights = cluster_weights,
                          standardized = standardized, method = method)
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

  # values for the from the previous check
  fn_old <- fn(w_new, abs_eps = abs_eps, rel_eps = rel_eps,
               minvls = minvls, maxvls = maxvls, indices = NULL)
  old_par <- w_new
  old_k <- k
  old_H <- H
  t_old <- t
  min_step_size <- step_factor / 2^8

  .check_if_increased <- function(){
    if(k %% check_every == 0){
      fn_new <- fn(w_new, abs_eps = abs_eps, rel_eps = rel_eps,
                   minvls = minvls, maxvls = maxvls, indices = NULL)
      var_est <- attr(fn_new, "std")^2 + attr(fn_old, "std")^2

      if(fn_new < fn_old + 2.326 * sqrt(var_est)){
        # fine as we are minimizing the negative likelihood
        fn_old <<- fn_new
        old_k <<- k
        old_H <<- H
        old_par <<- w_new
        t_old <<- t

      } else {
        # we reset
        H <<- old_H
        w_new <<- old_par
        step_factor <<- step_factor / 2
        t <<- t_old
        if(step_factor > min_step_size)
          # otherwise we can go on forever
          k <<- old_k - 1L

        if(verbose)
          cat(sprintf(
            "\nThe log-likelihood had not increased. The difference is %.6f from %.3f.\nHalving the step size to %f\n",
            -(fn_new - fn_old), -fn_old * n_terms, step_factor))

        return(FALSE)

      }
    }

    TRUE
  }

  while(k < n_it){
    w_vals[] <- NA_real_
    i_indices <- if(k < 1) 1L else 1:n_grad_steps
    for(i in i_indices){
      if((k <- k + 1L) > n_it)
        break

      # check if we need to reset
      if(!(is_ok <- .check_if_increased()))
        break

      # perform the gradient step
      w_old <- w_new
      S <- sample(indices, n_grad)
      gr_val <- gr(w_old, abs_eps = abs_eps, rel_eps = rel_eps,
                   minvls = minvls, maxvls = maxvls, indices = S)
      w_new <- w_old - step_factor / k * drop(H %*% gr_val)
      w_vals[, i] <- w_new
    }

    if(!is_ok)
      next

    if(i < min(n_grad_steps, length(i_indices)))
      # no need for an update of the Hessian approximation
      break

    t <- t + 1L
    omegas[, t] <- rowMeans(w_vals[, i_indices, drop = FALSE], na.rm = TRUE)
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
    s_y <- sum(s * y)
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

  omegas <- omegas[, apply(!is.na(omegas), 2L, all), drop = FALSE]
  list(par = w_new, omegas = omegas, H = H)
}


#' Compute Profile Likelihood Based Confidence Intervals
#' @inheritParams pedmod_opt
#'
#' @param par numeric vector with the maximum likelihood estimator e.g. from
#' \code{\link{pedmod_opt}}.
#' @param delta numeric scalar with an initial step to take. Subsequent steps
#' are taken by \code{2^(<iteration number> - 1) * delta}. Two times the
#' standard error is a good value or a guess thereof. Hessian approximations are
#' not implemented as of this writing and therefore the user needs to provide
#' some guess.
#' @param alpha numeric scalar with the confidence level required.
#' @param which_prof integer scalar with index of the parameter which the
#' profile likelihood curve should be computed for.
#' @param max_step integer scalar with the maximum number of steps to take in
#' either directions.
#' @param verbose logical for whether output should be printed to the console
#' during the estimation of the profile likelihood curve.
#' @param ... arguments passed on to \code{\link{pedmod_opt}}.
#'
#' @seealso
#' \code{\link{pedmod_opt}} and \code{\link{pedmod_sqn}}.
#'
#'
#' @return
#' A list with the following elements
#'   \item{confs}{2D numeric vector with the profile likelihood based confidence
#'                interval.}
#'   \item{xs}{the points at which the profile likelihood is evaluated.}
#'   \item{p_log_Lik}{the log profile likelihood values at \code{xs}.}
#'   \item{data}{list with the returned objects from \code{\link{pedmod_opt}}.}
#'
#' @importFrom stats spline approx qchisq qnorm setNames splinefun
#' @export
pedmod_profile <- function(ptr, par, delta, maxvls, minvls = -1L,
                           alpha = .05, abs_eps,
                           rel_eps, which_prof, indices = NULL,
                           do_reorder = TRUE, use_aprx = FALSE, n_threads = 1L,
                           cluster_weights = NULL, method = 0L, seed = 1L,
                           verbose = FALSE, max_step = 15L,
                           standardized = FALSE, ...){
  # checks
  stopifnot(
    is.numeric(par),
    length(which_prof) == 1L, is.integer(which_prof),
    which_prof %in% seq_along(par),
    is.numeric(delta), is.finite(delta), delta > 0,
    is.numeric(alpha), length(alpha) == 1, is.finite(alpha),
    alpha > 0, alpha < 1)

  # assign function to evaluate the log likelihood
  fn <- function(par, minv = minvls){
    set.seed(seed)
    eval_pedigree_ll(ptr = ptr, par = par, maxvls = maxvls, abs_eps = abs_eps,
                     rel_eps = rel_eps, indices = indices, minvls = minv,
                     do_reorder = do_reorder, use_aprx = use_aprx,
                     n_threads = n_threads, cluster_weights = cluster_weights,
                     standardized = standardized, method = method)
  }
  optim_res <- list(par = par, value = -fn(par))

  if(verbose)
    message(sprintf("The estimate of the standard error of the log likelihood is %.8f. Preferably this should be below 0.001",
                    attr(optim_res$value, "std")))

  # assign function to do the model fitting
  n_scales <- get_n_scales(ptr)
  total_var <- 1 + sum(exp(tail(par, n_scales)))
  beta_0 <- head(par, -n_scales) / sqrt(total_var)

  chi_val <- qchisq(1 - alpha, 1)
  crit_value <- -optim_res$value - chi_val / 2

  wrap_optim <- function(x, optim_obj, dir){
    if(verbose)
      message(sprintf("Log likelihood is %.4f at %f (critical value is %.4f)",
                      -optim_obj$value, x, crit_value))
    list(x = x, value = -optim_obj$value, optim = optim_obj,
         z_val = sign(dir) * sqrt((optim_obj$value - optim_res$value) * 2))
  }

  do_fit <- function(x, dir){
    # get the starting value
    par[which_prof] <- x
    if(which_prof > length(beta_0) && !standardized){
      total_var <- 1 + sum(exp(tail(par, n_scales)))
      par[seq_along(beta_0)] <- beta_0 * sqrt(total_var)
    }

    opt_out <- pedmod_opt(
      ptr = ptr, par = par, maxvls = maxvls, abs_eps = abs_eps, rel_eps = rel_eps,
      seed = seed, indices = indices, minvls = minvls, do_reorder = do_reorder,
      use_aprx = use_aprx, n_threads = n_threads, fix = which_prof,
      cluster_weights = cluster_weights, standardized = standardized,
      method = method, ...)

    wrap_optim(x, opt_out, dir)
  }

  # find points on the profile likelihood curve in either direction
  get_points <- function(dir){
    dir <- sign(dir)
    if(verbose)
      message(sprintf(
        "\nFinding the %s limit of the profile likelihood curve",
        if(dir < 0) "lower" else "upper"))

    step <- 0L
    out <- vector("list", max_step)
    prev <- -optim_res$value
    did_fail <- FALSE

    while(prev > crit_value && (step <- step + 1L) <= max_step){
      out[[step]] <- do_fit(par[which_prof] + dir *  2^(step - 1) * delta,
                            dir = dir)
      if(out[[step]]$value > prev){
        warning("Log likelihood did not decrease. Either 'optim_res' is not an optimum or the precision needs to be increased")
        did_fail <- TRUE
        break
      }

      prev <- out[[step]]$value
    }

    .report_failed <- function()
      if(verbose && step > max_step)
        warning(sprintf(
          "Failed to find the appropiate point in %d steps", max_step))

    .report_failed()

    if(did_fail || step > max_step)
      return(out[sapply(out, length) > 0])

    ub <- if(step == 1L)
      wrap_optim(par[which_prof], optim_res, 1) else out[[step - 1L]]
    lb <- out[[step]]

    while(ub$value - lb$value > chi_val / 6 && (step <- step + 1L) <= max_step){
      # compute the next value
      xs <- c(unlist(sapply(out, `[[`, "x"))    , par[which_prof])
      ys <- c(unlist(sapply(out, `[[`, "value")), -optim_res$value)
      sp <- splinefun(ys, xs, method = "monoH.FC")
      y <- seq(min(ys), max(ys), length.out = 4 * max_step)
      x <- sp(y)
      next_val <- approx(y, x, xout = crit_value)$y
      if(abs(next_val - ub$x) > abs(next_val - lb$x))
        next_val <- 8 / 9 * next_val + ub$x / 9
      else
        next_val <- 8 / 9 * next_val + lb$x / 9

      out[[step]] <- do_fit(next_val, dir = dir)

      if(out[[step]]$value > ub$value || out[[step]]$value < lb$value){
        warning("Log likelihood does not seem monotonic. Likely the precision needs to be increased")
        break
      }

      if(out[[step]]$value < crit_value)
        lb <- out[[step]]
      else
        ub <- out[[step]]
    }

    .report_failed()
    out <- out[sapply(out, length) > 0]
    out[order(sapply(out, `[[`, "x"))]
  }

  res_down <- get_points(-1)
  res_up   <- get_points( 1)
  out <- c(
    res_down, list(wrap_optim(par[which_prof], optim_res, dir = 0)), res_up)

  # compute the confidence interval
  xs  <- sapply(out, `[[`, "x")
  zs  <- sapply(out, `[[`, "z_val")
  pls <- sapply(out, `[[`, "value")
  sp <- spline(xs, zs)
  pvs <- c(alpha / 2, 1 - alpha/2)
  confs <- setNames(approx(sp$y, sp$x, xout = qnorm(pvs))$y,
                    sprintf("%.2f pct.", 100 * pvs))

  # return
  list(confs = confs, xs = xs, p_log_Lik = pls, data = out)
}
