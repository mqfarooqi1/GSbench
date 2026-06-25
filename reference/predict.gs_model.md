# Predict from a fitted genomic prediction model

Predict from a fitted genomic prediction model

## Usage

``` r
# S3 method for class 'gs_ensemble'
predict(object, newgeno, ...)

# S3 method for class 'gs_model'
predict(object, newgeno, ...)
```

## Arguments

- object:

  A `gs_model` from
  [`gs_fit()`](https://mqfarooqi1.github.io/GSbench/reference/gs_fit.md).

- newgeno:

  Marker matrix for the individuals to predict (same markers as
  training).

- ...:

  Ignored.

## Value

A numeric vector of predictions, one per row of `newgeno`.

## Examples

``` r
sim <- simulate_population(n = 120, m = 400, seed = 1)
fit <- gs_fit(sim$pheno[1:90], sim$geno[1:90, ], model = "gblup")
predict(fit, sim$geno[91:120, ])
#>  [1] -1.98866802 -0.29966079  1.78887571  2.48769719  0.04331597 -0.59946530
#>  [7]  2.21595315  3.58454879  1.51979193  1.77824768 -0.77063277  0.93959233
#> [13]  2.02117032 -1.74956335  3.02439391 -1.92682047  2.12224580  0.45306183
#> [19] -1.58324391 -0.81607869 -0.87959915 -0.28384698  1.97443102  0.78871045
#> [25]  0.30677387  2.01507714  2.48597472 -2.09485655 -0.89380108  1.74747468
```
