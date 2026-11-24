## Resubmission after archival

CRAN archived pedmod 0.2.4 on 2024-07-02. Version 0.2.5 fixes the invalid
references that GCC UBSan reported when graph functions used unit vertex
weights. It also fixes the Rd list for `pedigree_ll_terms()` and a Hessian
calculation that failed with current RcppArmadillo.

## Local checks

On Debian 13 with R 4.5.0, GCC 14.2.0 and RcppArmadillo 15.8.0-1,
`R CMD check --as-cran pedmod_0.2.5.tar.gz` completed with zero errors,
zero warnings and two notes. All 491 test assertions passed.

- The incoming feasibility check reports the previous archival.
- The check host could not verify the current time.

The installed package occupies 5.8 MB on this system. The two prebuilt
vignettes account for 3.9 MB, and the compiled library accounts for 1.6 MB.
This check reported the size as informational.
