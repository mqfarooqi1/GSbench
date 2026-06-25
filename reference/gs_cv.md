# Cross-validate genomic prediction models

Runs a breeding-relevant cross-validation and reports predictive ability
(the correlation between predictions and observed phenotypes on held-out
individuals) for each model. Two schemes are supported: `"kfold"`
(random k-fold, i.e. predicting untested lines, repeated `reps` times)
and `"leave_group_out"` (hold out one family/environment at a time, via
`groups`).

## Usage

``` r
gs_cv(
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

An object of class `gs_cv`: a list with `accuracy` (a data frame of
`model`, `mean`, `sd`, `n_folds`), `per_fold` (the raw fold results),
and the call settings.

## Details

Fitting is wrapped so that a model which errors on a fold records `NA`
for that fold rather than aborting the whole run.

## Examples

``` r
sim <- simulate_population(n = 120, m = 400, seed = 1)
cv <- gs_cv(sim$pheno, sim$geno, models = "gblup", k = 5, seed = 1)
cv
#> <gs_cv: kfold>
#>   5-fold x 1 rep(s)
#>  model  mean    sd n_folds
#>  gblup 0.153 0.106       5
#>   (accuracy = predictive ability, cor(pred, observed) on held-out data)
```
