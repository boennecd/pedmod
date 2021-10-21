#' @rdname eval_pedigree
#'
#' @title Approximate the Log Marginal Likelihood
#' @description
#' Approximate the log marginal likelihood and the derivatives with
#' respect to the model parameters.
#'
#' @param maxvls maximum number of samples in the approximation for each
#' marginal likelihood term.
#' @param minvls minimum number of samples for each
#' marginal likelihood term. Negative values provides a
#' default which depends on the dimension of the integration.
#' @param ptr object from \code{\link{pedigree_ll_terms}}.
#' @param par numeric vector with fixed effect coefficients and log scale
#' parameters. The log scale parameters should be last.
#' @param indices zero-based vector with indices of which log marginal
#' likelihood terms to include. Use \code{NULL} if all indices should be
#' used.
#' @param n_threads number of threads to use.
#' @param cluster_weights numeric vector with weights for each cluster. Use
#' \code{NULL} if all clusters have weight one.
#' @param standardized logical for whether to use the standardized or direct
#' parameterization. See \code{\link{standardized_to_direct}} and the vignette
#' at \code{vignette("pedmod", package = "pedmod")}.
#'
#' @inheritParams pedigree_ll_terms
#' @inheritParams mvndst
#'
#' @return \code{eval_pedigree_ll}:
#' a scalar with the log marginal likelihood approximation.
#' It has an attribute called \code{"n_fails"} which shows the number of
#' log marginal likelihood term approximations which do not satisfy
#' the \code{abs_eps} and \code{rel_eps} criteria and an attribute called
#' \code{std} with a standard error estimate based on the delta rule.
#'
#' @examples
#' # three families as an example
#' fam_dat <- list(
#'   list(
#'     y = c(FALSE, TRUE, FALSE, FALSE),
#'     X = structure(c(
#'       1, 1, 1, 1, 1.2922654151273, 0.358134905909256, -0.734963997107464,
#'       0.855235473516044, -1.16189500386223, -0.387298334620742,
#'       0.387298334620742, 1.16189500386223),
#'       .Dim = 4:3, .Dimnames = list( NULL, c("(Intercept)", "X1", ""))),
#'     rel_mat = structure(c(
#'       1, 0.5, 0.5, 0.125, 0.5, 1, 0.5, 0.125, 0.5, 0.5,
#'       1, 0.125, 0.125, 0.125, 0.125, 1), .Dim = c(4L, 4L)),
#'     met_mat = structure(c(1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1),
#'                         .Dim = c(4L, 4L))),
#'   list(
#'     y = c(FALSE, FALSE, FALSE),
#'     X = structure(c(
#'       1, 1, 1, -0.0388728997202442, -0.0913782435233639,
#'       -0.0801619722392612, -1, 0, 1), .Dim = c(3L, 3L)),
#'     rel_mat = structure(c(
#'       1, 0.5, 0.125, 0.5, 1, 0.125, 0.125, 0.125, 1), .Dim = c(3L, 3L)),
#'     met_mat = structure(c(
#'       1, 1, 0, 1, 1, 0, 0, 0, 1), .Dim = c(3L, 3L))),
#'   list(
#'     y = c(TRUE, FALSE),
#'     X = structure(c(
#'       1, 1, 0.305275750370738, -1.49482995913648,  -0.707106781186547,
#'       0.707106781186547),
#'       .Dim = 2:3, .Dimnames = list( NULL, c("(Intercept)", "X1", ""))),
#'     rel_mat = structure(c(1, 0.5,  0.5, 1), .Dim = c(2L, 2L)),
#'     met_mat = structure(c(1, 1, 1, 1), .Dim = c(2L,  2L))))
#'
#' # get the data into the format needed for the package
#' dat_arg <- lapply(fam_dat, function(x){
#'   # we need the following for each family:
#'   #   y: the zero-one outcomes.
#'   #   X: the design matrix for the fixed effects.
#'   #   scale_mats: list with the scale matrices for each type of effect.
#'   list(y = as.numeric(x$y), X = x$X,
#'        scale_mats = list(x$rel_mat, x$met_mat))
#' })
#'
#' # get a pointer to the C++ object
#' ptr <- pedigree_ll_terms(dat_arg, max_threads = 1L)
#'
#' # approximate the log marginal likelihood
#' beta <- c(-1, 0.3, 0.2) # fixed effect coefficients
#' scs <- c(0.5, 0.33)     # scales parameters
#'
#' set.seed(44492929)
#' system.time(ll1 <- eval_pedigree_ll(
#'   ptr = ptr, par = c(beta, log(scs)), abs_eps = -1, maxvls = 1e5,
#'   rel_eps = 1e-5, minvls = 2000, use_aprx = FALSE))
#' ll1 # the approximation
#'
#' # with the approximation of pnorm and qnorm
#' system.time(ll2 <- eval_pedigree_ll(
#'   ptr = ptr, par = c(beta, log(scs)), abs_eps = -1, maxvls = 1e5,
#'   rel_eps = 1e-5, minvls = 2000, use_aprx = TRUE))
#' all.equal(ll1, ll2, tolerance = 1e-5)
#'
#' # cluster weights can be used as follows to repeat the second family three
#' # times and remove the third
#' system.time(deriv_w_weight <- eval_pedigree_grad(
#'   ptr = ptr, par = c(beta, log(scs)), abs_eps = -1, maxvls = 1e6,
#'   rel_eps = 1e-3, minvls = 2000, use_aprx = TRUE,
#'   cluster_weights = c(1, 3, 0)))
#'
#' # the same as manually repeating second cluster and not including the third
#' dum_dat <- dat_arg[c(1, 2, 2, 2)]
#' dum_ptr <- pedigree_ll_terms(dum_dat, 1L)
#' system.time(deriv_dum <- eval_pedigree_grad(
#'   ptr = dum_ptr, par = c(beta, log(scs)), abs_eps = -1, maxvls = 1e6,
#'   rel_eps = 1e-3, minvls = 2000, use_aprx = TRUE))
#' all.equal(deriv_dum, deriv_w_weight, tolerance = 1e-3)
#'
#' @export
eval_pedigree_ll <- function(ptr, par, maxvls, abs_eps, rel_eps,
                             indices = NULL, minvls = -1L, do_reorder = TRUE,
                             use_aprx = FALSE, n_threads = 1L,
                             cluster_weights = NULL, standardized = FALSE,
                             method = 0L){
  if(standardized)
    par <- standardized_to_direct(par = par, n_scales = get_n_scales(ptr))
  eval_pedigree_ll_cpp(
    ptr = ptr, par = par, maxvls = maxvls, abs_eps = abs_eps, rel_eps = rel_eps,
    indices = indices, minvls = minvls, do_reorder = do_reorder,
    use_aprx = use_aprx, n_threads = n_threads,
    cluster_weights = cluster_weights, method = method)
}

#' @rdname eval_pedigree
#'
#' @return \code{eval_pedigree_grad}: a vector with the derivatives with
#' respect to \code{par}. An attribute called \code{"logLik"} contains the
#' log marginal likelihood approximation. There will also be \code{"n_fails"}
#' attribute like for \code{eval_pedigree_ll} and an attribute called
#' \code{"std"} which first element is the standard error estimate of the
#' log likelihood based on the delta method and the last elements are the
#' standard error estimates of the gradient. The latter ignores the Monte Carlo
#' error from the likelihood approximation.
#'
#' @export
eval_pedigree_grad <- function(ptr, par, maxvls, abs_eps, rel_eps,
                               indices = NULL, minvls = -1L, do_reorder = TRUE,
                               use_aprx = FALSE, n_threads = 1L,
                               cluster_weights = NULL, standardized = FALSE,
                               method = 0L){
  if(standardized)
    par <- standardized_to_direct(par = par, n_scales = get_n_scales(ptr),
                                  jacobian = TRUE)

  out <- eval_pedigree_grad_cpp(
    ptr = ptr, par = par, maxvls = maxvls, abs_eps = abs_eps, rel_eps = rel_eps,
    indices = indices, minvls = minvls, do_reorder = do_reorder,
    use_aprx = use_aprx, n_threads = n_threads,
    cluster_weights = cluster_weights, method = method)

  if(!standardized)
    return(out)

  res <- drop(out %*% attr(par, "jacobian"))
  attributes(res) <- attributes(out)
  res
}
