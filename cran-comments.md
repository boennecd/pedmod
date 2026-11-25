## New version submission

CRAN archived pedmod 0.2.4 on 2024-07-02. Version 0.2.5 fixes the invalid
references that GCC UBSan reported when graph functions used unit vertex
weights.

## Check environments

- Debian 13 x86_64, R 4.5.0, GCC 14.2.0, `R CMD check --as-cran`.
- Ubuntu 26.04.1 x86_64, R 4.5.3, 4.6.2, and R-devel, GCC 15.2.0.
  macOS 26.7.2 arm64, R 4.6.2, Apple Clang 21.0.0.
  Windows Server 2022 x64, R 4.6.2, GCC 14.3.0.
  All used `R CMD check --as-cran --no-manual`.
- Fedora 44 x86_64, R-devel, GCC 16.2.1 with
  `-fsanitize=address,undefined,bounds-strict` and
  `R CMD check --as-cran --no-manual --no-build-vignettes`.
- Ubuntu 26.04.1 x86_64, R-devel, Clang 23.1.6 with ASan and UBSan.
  Both checks used `R CMD check --no-manual --no-stop-on-test-error`.
- Fedora 44 x86_64, R-devel, GCC 16.2.1, Valgrind 3.27.1,
  `R CMD check --no-manual --use-valgrind --no-stop-on-test-error`.
- macOS Sequoia 1.52 with a M1 chip.
- devtools::check_win_devel(), devtools::check_win_release() and 
  devtools::check_win_oldrelease().

## Results

The standard checks reported no errors or warnings. When there is a note then it 
is about the prior archival.
