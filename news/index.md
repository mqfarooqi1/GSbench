# Changelog

## GSbench 0.1.0

CRAN release: 2026-06-30

First release.

- Core (base R):
  [`simulate_population()`](https://mqfarooqi1.github.io/GSbench/reference/simulate_population.md),
  [`qc_markers()`](https://mqfarooqi1.github.io/GSbench/reference/qc_markers.md),
  [`impute_markers()`](https://mqfarooqi1.github.io/GSbench/reference/impute_markers.md),
  [`Gmatrix()`](https://mqfarooqi1.github.io/GSbench/reference/Gmatrix.md)
  (VanRaden), and
  [`gblup()`](https://mqfarooqi1.github.io/GSbench/reference/gblup.md)
  (GBLUP by REML, validated against
  [`rrBLUP::mixed.solve`](https://rdrr.io/pkg/rrBLUP/man/mixed.solve.html)).
- Unified modelling interface
  [`gs_fit()`](https://mqfarooqi1.github.io/GSbench/reference/gs_fit.md)
  / [`predict()`](https://rdrr.io/r/stats/predict.html) covering GBLUP,
  elastic net (`glmnet`), random forest (`ranger`) and gradient boosting
  (`xgboost`).
- [`gs_cv()`](https://mqfarooqi1.github.io/GSbench/reference/gs_cv.md)
  for breeding-relevant cross-validation (random k-fold and
  leave-one-group-out).
- [`gs_ensemble()`](https://mqfarooqi1.github.io/GSbench/reference/gs_ensemble.md),
  a stacked super-learner combining base models with non-negative,
  out-of-fold-fitted weights.
- [`gs_benchmark()`](https://mqfarooqi1.github.io/GSbench/reference/gs_benchmark.md)
  with `print`/`summary`/`plot` to compare all available models under
  one cross-validation.
