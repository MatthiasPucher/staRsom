## R CMD check results

Local check (`R CMD check --as-cran --run-donttest`, CRAN incoming checks without remote access):

0 errors | 0 warnings | 0 notes

* This is a new submission.
* Words in DESCRIPTION that a spell check may report are abbreviations (DOM, EEM, EEMs, PARAFAC, SOM) and names (Kohonen, Ejarque, Butturini).

## Test environments

* local: Ubuntu 24.04, R 4.6.1, single-threaded BLAS
* GitHub Actions: (to be added after the first run)
* win-builder, R-devel: (to be added)

## Notes for the reviewers

* Random numbers: functions using random numbers have an argument `seed`. The default `NULL` uses the current random number stream of the session. Where seeds are set internally (a number given as `seed`, and the seed of each map in `eemsom_stability()`), the previous state of the random number generator is restored afterwards by assigning `.Random.seed` in the global environment, like `withr::with_seed()` does. Nothing else in the global environment is changed.
* Parallel computing: examples and tests use at most 2 cores.
* With a multi-threaded OpenBLAS, the CPU time of the tests and of some examples exceeds 2.5 times the elapsed time. This is caused by the BLAS library (matrix products and singular value decompositions), not by the package; with a single BLAS thread, no example and no test exceeds the limits for CPU time and elapsed time.
