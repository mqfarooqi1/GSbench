# Predict breeding values for new genotypes from a GBLUP fit

Predict breeding values for new genotypes from a GBLUP fit

## Usage

``` r
# S3 method for class 'gblup'
predict(object, newgeno, ...)
```

## Arguments

- object:

  A `gblup` fit created with `geno` supplied.

- newgeno:

  Marker matrix for new individuals (markers must match the training
  markers, in the same order).

- ...:

  Ignored.

## Value

A numeric vector of predicted values (intercept + marker effects), one
per row of `newgeno`.

## Examples

``` r
sim <- simulate_population(n = 150, m = 400, seed = 1)
fit <- gblup(sim$pheno[1:120], sim$geno[1:120, ])
predict(fit, sim$geno[121:150, ])
#>  [1]  2.2394786 -0.1549250  1.5053554  2.5441656  4.5218152 -1.7236605
#>  [7]  2.5763509  2.1231666  1.7740676  0.4868161 -0.7821772  2.7136004
#> [13]  1.9608765  2.0175759 -4.6569745  2.5057748 -0.5546527 -3.5618100
#> [19]  0.5803449  2.3974825 -4.1287929  3.4608844  5.8706409  1.2929954
#> [25]  4.5799813  2.8146889  0.3444057  0.2098574  2.9769580  3.6270088
```
