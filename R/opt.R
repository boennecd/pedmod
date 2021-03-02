#' Optimize the Log Marginal Likelihood
#'
#' Optimizes \code{\link{eval_pedigree_ll}} and \code{\link{eval_pedigree_grad}}
#' using a passed optimization function.
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
#' The output from the \code{opt_func} argument.
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
