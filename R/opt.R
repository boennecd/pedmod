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
#' for computing profile likelihoods. \code{NULL} yield all parameters.
#' @param ... Arguments passed to \code{opt_func}.
#'
#' @return
#' \code{pedmod_opt}: The output from the \code{opt_func} argument.
#'
#' @importFrom stats optim
#' @export
pedmod_opt <- function(ptr, par, maxvls, abs_eps, rel_eps,
                       opt_func = NULL, seed = 1L, indices = NULL, minvls = -1L,
                       do_reorder = TRUE, use_aprx = FALSE, n_threads = 1L,
                       cluster_weights = NULL, fix = NULL, ...){
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
      cluster_weights = cluster_weights), silent = TRUE)
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
                          cluster_weights = cluster_weights)
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
                         opt_func = NULL, seed = 1L, indices = NULL,
                         minvls = 100L, do_reorder = TRUE, use_aprx = TRUE,
                         n_threads = 1L, cluster_weights = NULL, ...){
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
  logLik_no_rng <- -sum(start_fit$deviance) / 2
  n_scales <- length(data[[1L]]$scale_mats)
  succes <- FALSE
  for(i in 1:10){
    # there seems to be a tendency to overstep. This __may__ solve it
    sc <- rep(log(.1 * i^2), n_scales)
    stopifnot(length(sc) > 0)
    opt <- optim(sc, fn, gr, method = "BFGS")

    if(sum(exp(opt$par)) < 100){
      # found a solution where the total variance is not -> infity
      sc <- opt$par
      beta_scaled <- beta * sqrt(1 + sum(exp(sc)))
      succes <- TRUE
      logLik_est <- -opt$value
      break
    }
  }

  if(!succes){
    # fall back to use low values for the scale parameters as the starting
    # values
    sc <- rep(log(.01), n_scales)
    beta_scaled <- beta * sqrt(1 + sum(exp(sc)))
    logLik_est <- -fn(sc)
  }

  # return
  list(par = c(beta_scaled, sc), beta_no_rng = beta,
       logLik_no_rng = logLik_no_rng,
       logLik_est = -opt$value)
}
