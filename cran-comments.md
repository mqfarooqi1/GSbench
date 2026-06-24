## Submission

New submission: GSbench 0.1.0 — a unified interface to fit, cross-validate and
benchmark genomic-selection models (GBLUP and ridge marker effects in base R,
plus optional machine-learning backends and a stacked ensemble) from SNP data.

## R CMD check results

Local `R CMD check --as-cran` (Windows 11, R 4.5.2): 0 errors | 0 warnings |
1 note (the standard "New submission").

## Notes

* The GBLUP solver is implemented in base R and validated against
  `rrBLUP::mixed.solve` in the test suite (identical variance components; GEBVs
  correlating at ~1).
* Machine-learning backends (`glmnet`, `ranger`, `xgboost`) are optional
  (Suggests) and used only when installed; tests and examples skip gracefully if
  they are absent.
* Functions do not write to the user's file system or modify the global
  environment; the `seed` arguments use `withr::with_seed()`.
* Method references are given in the Description with DOIs.

## Test environments

* Local: Windows 11, R 4.5.2
* win-builder (release + devel) planned before release
