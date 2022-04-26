# pedmod 0.2.0
* `pedigree_ll_terms_loadings` is implemented to support models with individual 
  specific covariance scale parameters (e.g. individual specific 
  heritabilities).
* The minimax tilting method suggested by Botev (2017) (see 
  https://doi.org/10.1111/rssb.12162) is implemented. The method is less 
  numerically stable and thus required more care when implementing. This yield a 
  higher per randomized quasi-Monte Carlo sample cost. Though, the increased 
  cost may be worthwhile for low probability events because of a reduced 
  variance at a fixed number of samples.
* The `vls_scales` argument is added which allows the user to use more 
  randomized quasi-Monte Carlo samples for some log likelihood terms. This is 
  useful e.g. when one uses weighted terms.

# pedmod 0.1.0 
* First release on CRAN.
