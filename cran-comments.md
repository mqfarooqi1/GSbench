## Resubmission

This resubmission addresses the NOTEs from the automated incoming pre-test.

### "Possibly misspelled words in DESCRIPTION"
These are not misspellings:

* BLUP, GBLUP - established acronyms in quantitative genetics (genomic best
  linear unbiased prediction), already spelled out in the Description.
* Endelman, Meuwissen, VanRaden - author surnames in the method references.
* Benchmarking - a correct English word (the first word of the Title).

No change was needed for these; they are flagged because they are domain
acronyms and proper names.

### "CPU time N times elapsed time" in tests (Debian)
Fixed. The affected tests exercised the optional machine-learning backends
(ranger, xgboost), which use OpenMP threads internally. Those tests are now
guarded with `testthat::skip_on_cran()`, so CRAN's test run stays single-core;
they continue to run locally. The remaining on-CRAN tests use only base R and
the GBLUP solver.

## R CMD check results

Local `R CMD check --as-cran` (Windows 11, R 4.5.2): 0 errors | 0 warnings |
1 note (the standard "New submission").

## Notes

* The GBLUP solver is base R and validated against `rrBLUP::mixed.solve` in the
  (non-CRAN) tests.
* Machine-learning backends are optional (Suggests) and thread-limited
  (`num.threads = 1`, `nthread = 1`); tests using them now skip on CRAN.
* Functions write only to `tempdir()` and do not modify the global environment
  (seeds via `withr::with_seed()`). Method references carry DOIs.

## Test environments

* Local: Windows 11, R 4.5.2
* win-builder (r-devel): the pre-test NOTEs above, now addressed.
