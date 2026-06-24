# GSbench 0.1.0

First release.

* Core (base R): `simulate_population()`, `qc_markers()`, `impute_markers()`,
  `Gmatrix()` (VanRaden), and `gblup()` (GBLUP by REML, validated against
  `rrBLUP::mixed.solve`).
* Unified modelling interface `gs_fit()` / `predict()` covering GBLUP, elastic
  net (`glmnet`), random forest (`ranger`) and gradient boosting (`xgboost`).
* `gs_cv()` for breeding-relevant cross-validation (random k-fold and
  leave-one-group-out).
* `gs_ensemble()`, a stacked super-learner combining base models with
  non-negative, out-of-fold-fitted weights.
* `gs_benchmark()` with `print`/`summary`/`plot` to compare all available
  models under one cross-validation.
