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

## Remote checks

GitHub Actions passed all five standard checks on commit `175e4dd`.
The checks covered R release on Linux, macOS, and Windows, plus R-devel and
R oldrel-1 on Linux. The R-devel check reported no errors, warnings, or notes.
[GitHub Actions run](https://github.com/boennecd/pedmod/actions/runs/45491678393)

On commit `1092da6`, R-hub's `gcc-asan` check ran with GCC AddressSanitizer
and UndefinedBehaviorSanitizer. The examples and tests passed without the
reported misaligned-reference diagnostic or other sanitizer findings.
The check reported one note about deprecated matrix attribute names.
Commit `175e4dd` updates those names in the examples and tests.
The `clang-asan` and `clang-ubsan` checks also passed.
[R-hub sanitizer run](https://github.com/boennecd/pedmod/actions/runs/45489278065)

On the same commit, `R CMD check` reported status `OK` under Valgrind.
R-hub marked the job as failed after Valgrind reported one 384-byte
allocation as possibly lost in GNU OpenMP worker-thread setup. Valgrind
reported zero bytes definitely or indirectly lost. A standalone
two-thread OpenMP program reproduced the same allocation stack without
loading pedmod.
[Valgrind log](https://github.com/boennecd/pedmod/actions/runs/45489278065/job/138783022246)
