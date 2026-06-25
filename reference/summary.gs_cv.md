# Summary of a cross-validation / benchmark

Summary of a cross-validation / benchmark

## Usage

``` r
# S3 method for class 'gs_cv'
summary(object, ...)
```

## Arguments

- object:

  A `gs_cv` object.

- ...:

  Ignored.

## Value

The accuracy data frame, invisibly.

## Examples

``` r
sim <- simulate_population(n = 100, m = 300, seed = 1)
summary(gs_cv(sim$pheno, sim$geno, models = "gblup", seed = 1))
#> <gs_cv: kfold>
#>   5-fold x 1 rep(s)
#>  model mean    sd n_folds
#>  gblup 0.11 0.276       5
#>   (accuracy = predictive ability, cor(pred, observed) on held-out data)
```
