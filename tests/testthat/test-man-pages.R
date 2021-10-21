context("testing examples in manual pages")

test_that("examples in manual pages gives the correct answer for eval_pedigree_[ll]/[grad]", {
  S <- structure(c(0.766725077267174, 0.249137157109195, -0.301516767233159, 0.112735950556606, -0.121636706660597, 0.0853052478227257, -0.0111813042831139, -0.0201254731753443, 0.149427110129624, -0.021996990666367, 0.249137157109195, 1.42410274977015, -0.338609740421468, -0.0425085949907385, -0.613458260081048, -0.292986225195268, 0.126875549501989, 0.0939246378491726, 0.00591847504584558, 0.221189267455151, -0.301516767233159, -0.338609740421468, 1.13687004435125, 0.303669017834393, 0.399060537444824, -0.0256113853532709, 0.220769628380184, -0.0219655194790159, -0.107070143975143, 0.0556497760193132, 0.112735950556606, -0.0425085949907385, 0.303669017834393, 0.854325399612162, 0.0924461969299186, 0.0860338995123319, -0.058252637830639, -0.274008787179135, 0.137462284747638, -0.102623896963213, -0.121636706660597, -0.613458260081048, 0.399060537444824, 0.0924461969299186, 1.18847861341054, 0.137701677167591, -0.257837324154815, -0.125751791453519, 0.0941660348081916, 0.000446561075573518, 0.0853052478227257, -0.292986225195268, -0.0256113853532709, 0.0860338995123319, 0.137701677167591, 0.792314183442372, -0.0574683130409049, -0.127914440182539, -0.0819510156005461, -0.281344369049677, -0.0111813042831139, 0.126875549501989, 0.220769628380184, -0.058252637830639, -0.257837324154815, -0.0574683130409049, 1.0494755855744, 0.0680630636010323, -0.22917987710939, 0.339801755016208, -0.0201254731753443, 0.0939246378491726, -0.0219655194790159, -0.274008787179135, -0.125751791453519, -0.127914440182539, 0.0680630636010323, 0.66321514703524, 0.00981756084965376, 0.396173748919282, 0.149427110129624, 0.00591847504584558, -0.107070143975143, 0.137462284747638, 0.0941660348081916, -0.0819510156005461, -0.22917987710939, 0.00981756084965376, 0.967089365130523, -0.0792898297089128, -0.021996990666367, 0.221189267455151, 0.0556497760193132, -0.102623896963213, 0.000446561075573518, -0.281344369049677, 0.339801755016208, 0.396173748919282, -0.0792898297089128, 1.12782790972313), .Dim = c(10L, 10L))
  u <- c(-1.04413462631653, 0.569719627442413, -0.135054603880824, 2.40161776050478, -0.0392400027331692, 0.689739362450777, 0.0280021587806661, -0.743273208882405, 0.188792299514343, -1.80495862889104)

  n <- NCOL(S)
  # truth <- mvndst(
  #   lower = rep(-Inf, n), upper = u, sigma = S, mu = numeric(n),
  #   maxvls = 1e8, abs_eps = 0, rel_eps = 1e-6, use_aprx = FALSE)
  truth <- structure(
    7.34188854727007e-05, n_it = 8514544L, inform = 0L,
    abserr = 6.991581314015e-11)

  set.seed(1)
  pedmod_res <- mvndst(
    lower = rep(-Inf, n), upper = u, sigma = S, mu = numeric(n),
    maxvls = 1e6, abs_eps = 0, rel_eps = 1e-4, use_aprx = TRUE)
  expect_equal(truth, pedmod_res, check.attributes = FALSE,
               tolerance = 1e-4)

  pedmod_res <- mvndst(
    lower = rep(-Inf, n), upper = u, sigma = S, mu = numeric(n),
    maxvls = 1e6, abs_eps = 0, rel_eps = 1e-4, use_aprx = TRUE,
    method = 1L)
  expect_equal(truth, pedmod_res, check.attributes = FALSE,
               tolerance = 1e-4)
})

test_that("examples in manual pages gives the correct answer for eval_pedigree_[ll]/[grad]", {
  # three families as an example
  fam_dat <- list(
    list(
      y = c(FALSE, TRUE, FALSE, FALSE),
      X = structure(c(1, 1, 1, 1, 1.2922654151273, 0.358134905909256, -0.734963997107464, 0.855235473516044, -1.16189500386223, -0.387298334620742, 0.387298334620742, 1.16189500386223), .Dim = 4:3, .Dimnames = list( NULL, c("(Intercept)", "X1", ""))),
      rel_mat = structure(c(1, 0.5, 0.5, 0.125, 0.5, 1, 0.5, 0.125, 0.5, 0.5, 1, 0.125, 0.125, 0.125, 0.125, 1), .Dim = c(4L, 4L)),
      met_mat = structure(c(1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1), .Dim = c(4L, 4L))),
    list(
      y = c(FALSE, FALSE, FALSE),
      X = structure(c(1, 1, 1, -0.0388728997202442, -0.0913782435233639, -0.0801619722392612, -1, 0, 1), .Dim = c(3L, 3L)),
      rel_mat = structure(c(1, 0.5, 0.125, 0.5, 1, 0.125, 0.125, 0.125, 1), .Dim = c(3L, 3L)),
      met_mat = structure(c(1, 1, 0, 1, 1, 0, 0, 0, 1), .Dim = c(3L, 3L))),
    list(
      y = c(TRUE, FALSE),
      X = structure(c(1, 1, 0.305275750370738, -1.49482995913648,  -0.707106781186547, 0.707106781186547), .Dim = 2:3, .Dimnames = list( NULL, c("(Intercept)", "X1", ""))),
      rel_mat = structure(c(1, 0.5,  0.5, 1), .Dim = c(2L, 2L)),
      met_mat = structure(c(1, 1, 1, 1), .Dim = c(2L,  2L))))

  # get the data into the format needed for the package
  dat_arg <- lapply(fam_dat, function(x){
    # we need the following for each family:
    #   y: the zero-one outcomes.
    #   X: the design matrix for the fixed effects.
    #   scale_mats: list with the scale matrices for each type of effect.
    list(y = as.numeric(x$y), X = x$X,
         scale_mats = list(x$rel_mat, x$met_mat))
  })

  # get a pointer to the C++ object
  ptr <- pedigree_ll_terms(dat_arg, max_threads = 1L)

  # approximate the log marginal likelihood
  beta <- c(-1, 0.3, 0.2) # fixed effect coefficients
  scs <- c(0.5, 0.33)

  # truth <- eval_pedigree_ll(
  #   ptr = ptr, par = c(beta, log(scs)), abs_eps = -1, maxvls = 1e9,
  #   rel_eps = 1e-8, minvls = 2000, use_aprx = FALSE)
  truth <- structure(-5.30140009701486, n_fails = 0L,
                     std = 3.3296139122967e-09)

  set.seed(44492929)
  ll1 <- eval_pedigree_ll(
    ptr = ptr, par = c(beta, log(scs)), abs_eps = -1, maxvls = 1e6,
    rel_eps = 1e-5, minvls = 2000, use_aprx = FALSE)
  expect_equal(ll1, truth, tolerance = 1e-5)

  # with the approximation of pnorm and qnorm
  ll2 <- eval_pedigree_ll(
    ptr = ptr, par = c(beta, log(scs)), abs_eps = -1, maxvls = 1e6,
    rel_eps = 1e-5, minvls = 2000, use_aprx = TRUE)
  expect_equal(ll2, truth, tolerance = 1e-5)

  # with Sobol sequences
  ll3 <- eval_pedigree_ll(
    ptr = ptr, par = c(beta, log(scs)), abs_eps = -1, maxvls = 1e6,
    rel_eps = 1e-5, minvls = 2000, use_aprx = FALSE, method = 1L)
  expect_equal(ll3, truth, tolerance = 1e-5)

  # w/ weights
  # truth_dat <- dat_arg[c(1, 2, 2, 2)]
  # truth_ptr <- pedigree_ll_terms(truth_dat, 1L)
  # deriv_truth <- eval_pedigree_grad(
  #   ptr = truth_ptr, par = c(beta, log(scs)), abs_eps = -1, maxvls = 1e8,
  #   rel_eps = 1e-6, minvls = 2000, use_aprx = FALSE)
  deriv_truth <- structure(
    c(-2.39509630042317, -0.10824194375542, -0.940601039742817,
      -0.314925453459061, -0.278867316602556),
    logLik = -4.88592585105156, n_fails = 0L,
    std = c(1.17013082330889e-09, 4.84625759211153e-08, 2.05599876859111e-08,
            3.71573106370985e-08, 9.95606406437179e-08, 7.01965244943291e-08))

  deriv_w_weight <- eval_pedigree_grad(
    ptr = ptr, par = c(beta, log(scs)), abs_eps = -1, maxvls = 1e6,
    rel_eps = 1e-3, minvls = 2000, use_aprx = TRUE,
    cluster_weights = c(1, 3, 0))
  expect_equal(deriv_w_weight, deriv_truth, tolerance = 1e-3)

  ll_w_weight <- eval_pedigree_ll(
    ptr = ptr, par = c(beta, log(scs)), abs_eps = -1, maxvls = 1e6,
    rel_eps = 1e-3, minvls = 2000, use_aprx = TRUE,
    cluster_weights = c(1, 3, 0))
  expect_equal(c(ll_w_weight), attr(deriv_truth, "logLik"), tolerance = 1e-3)
})
