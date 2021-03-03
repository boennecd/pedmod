% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opt.R
\name{pedmod_opt}
\alias{pedmod_opt}
\alias{pedmod_start}
\title{Optimize the Log Marginal Likelihood}
\usage{
pedmod_opt(
  ptr,
  par,
  maxvls,
  abs_eps,
  rel_eps,
  opt_func = NULL,
  seed = 1L,
  indices = NULL,
  minvls = -1L,
  do_reorder = TRUE,
  use_aprx = FALSE,
  n_threads = 1L,
  cluster_weights = NULL,
  fix = NULL,
  ...
)

pedmod_start(
  ptr,
  data,
  maxvls = 1000L,
  abs_eps = 0,
  rel_eps = 0.01,
  opt_func = NULL,
  seed = 1L,
  indices = NULL,
  minvls = 100L,
  do_reorder = TRUE,
  use_aprx = TRUE,
  n_threads = 1L,
  cluster_weights = NULL,
  ...
)
}
\arguments{
\item{ptr}{object from \code{\link{get_pedigree_ll_terms}}.}

\item{par}{starting values passed to \code{opt_func}.}

\item{maxvls}{maximum number of samples in the approximation.}

\item{abs_eps}{absolute convergence threshold.}

\item{rel_eps}{relative convergence threshold.}

\item{opt_func}{function to perform minimization with arguments like
\code{\link{optim}}. BFGS is used with \code{\link{optim}} if this argument
is \code{NULL}.}

\item{seed}{seed to pass to \code{\link{set.seed}} before each gradient and
function evaluation. Use \code{NULL} if the seed should not be fixed.}

\item{indices}{zero-based vector with indices of which log marginal
likelihood terms to include. Use \code{NULL} if all indices should be
used.}

\item{minvls}{minimum number of samples. Negative values provides a
default which depends on the dimension of the integration.}

\item{do_reorder}{\code{TRUE} if a heuristic variable reordering should
be used. \code{TRUE} is likely the best value.}

\item{use_aprx}{\code{TRUE} if a less precise approximation of
\code{\link{pnorm}} and \code{\link{qnorm}} should be used. This may
reduce the computation time while not affecting the result much.}

\item{n_threads}{number of threads to use.}

\item{cluster_weights}{numeric vector with weights for each cluster. Use
\code{NULL} if all clusters have weight one.}

\item{fix}{integer vector with indices of \code{par} to fix. This is useful
for computing profile likelihoods. \code{NULL} yield all parameters.}

\item{...}{Arguments passed to \code{opt_func}.}

\item{data}{the \code{\link{list}} that was passed to
\code{\link{get_pedigree_ll_terms}}.}
}
\value{
\code{pedmod_opt}: The output from the \code{opt_func} argument.

\code{pedmod_start}: A \code{list} with:
\itemize{
 \item par: the starting value.
 \item beta_no_rng: the fixed effects MLEs without random effects.
 \item logLik_no_rng: the log maximum likelihood without random effects.
 \item logLik_est: the likelihood at par.
}
}
\description{
Optimizes \code{\link{eval_pedigree_ll}} and \code{\link{eval_pedigree_grad}}
using a passed optimization function.
}
\details{
\code{pedmod_start} yield starting values which can be used for
\code{pedmod_opt}. The method is based on a heuristic where we make the
assumption that the fixed effects are unrelated to the random effects.
}