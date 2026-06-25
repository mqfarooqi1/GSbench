# Fit a genomic prediction model

One interface for several model families: `"gblup"` (always available),
`"elastic_net"` (glmnet), `"random_forest"` (ranger), `"xgboost"`
(xgboost), and `"ensemble"` (a stacked super-learner; see
[`gs_ensemble()`](https://mqfarooqi1.github.io/GSbench/reference/gs_ensemble.md)).
The returned object has a
[predict()](https://mqfarooqi1.github.io/GSbench/reference/predict.gs_model.md)
method that takes a new marker matrix.

## Usage

``` r
gs_fit(y, geno, model = "gblup", ...)
```

## Arguments

- y:

  Numeric phenotype vector (length n), no missing values.

- geno:

  Marker matrix (n x m, coded 0/1/2, no missing values).

- model:

  Model name; see
  [`available_models()`](https://mqfarooqi1.github.io/GSbench/reference/available_models.md).

- ...:

  Model-specific hyperparameters (e.g. `alpha` for elastic net,
  `num.trees` for random forest, `nrounds`/`eta`/`max_depth` for
  xgboost, `base_models` for the ensemble).

## Value

An object of class `gs_model` wrapping the fitted model.

## Examples

``` r
sim <- simulate_population(n = 120, m = 400, seed = 1)
fit <- gs_fit(sim$pheno, sim$geno, model = "gblup")
head(predict(fit, sim$geno))
#> [1] -3.8752143 -5.4217031  0.7654262  6.0236117  5.7933099 -4.5519275
```
