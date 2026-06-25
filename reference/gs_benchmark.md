# Benchmark all available genomic prediction models

Convenience wrapper around
[`gs_cv()`](https://mqfarooqi1.github.io/GSbench/reference/gs_cv.md)
that evaluates every model available in the session (including the
stacked `"ensemble"`) under one cross-validation, so they can be
compared on an equal footing. Returns a `gs_cv` object; use
[plot()](https://mqfarooqi1.github.io/GSbench/reference/plot.gs_cv.md)
to draw the comparison.

## Usage

``` r
gs_benchmark(
  y,
  geno,
  models = available_models(),
  k = 5,
  reps = 1,
  scheme = c("kfold", "leave_group_out"),
  groups = NULL,
  seed = NULL,
  ...
)
```

## Arguments

- y:

  Numeric phenotype vector (length n), no missing values.

- geno:

  Marker matrix (n x m, 0/1/2, no missing values).

- models:

  Models to evaluate; defaults to all available (see
  [`available_models()`](https://mqfarooqi1.github.io/GSbench/reference/available_models.md)).
  Unavailable models are dropped with a warning.

- k:

  Number of folds for `"kfold"`. Default 5.

- reps:

  Number of repeats for `"kfold"`. Default 1.

- scheme:

  `"kfold"` or `"leave_group_out"`.

- groups:

  Grouping vector (length n) for `"leave_group_out"`.

- seed:

  Optional seed (applied via
  [`withr::with_seed()`](https://withr.r-lib.org/reference/with_seed.html)).

- ...:

  Passed to
  [`gs_fit()`](https://mqfarooqi1.github.io/GSbench/reference/gs_fit.md)
  (model hyperparameters).

## Value

A `gs_cv` object (see
[`gs_cv()`](https://mqfarooqi1.github.io/GSbench/reference/gs_cv.md)).

## Examples

``` r
sim <- simulate_population(n = 120, m = 400, seed = 1)
bench <- gs_benchmark(sim$pheno, sim$geno, models = "gblup", k = 5, seed = 1)
bench$accuracy
#>   model      mean        sd n_folds
#> 1 gblup 0.1530038 0.1057235       5
```
