# Plot a model comparison

Bar chart of mean predictive ability per model, with +/- 1 SD whiskers
across folds.

## Usage

``` r
# S3 method for class 'gs_cv'
plot(x, ...)
```

## Arguments

- x:

  A `gs_cv` object.

- ...:

  Passed to
  [`graphics::barplot()`](https://rdrr.io/r/graphics/barplot.html).

## Value

`x`, invisibly.

## Examples

``` r
sim <- simulate_population(n = 120, m = 400, seed = 1)
# \donttest{
plot(gs_cv(sim$pheno, sim$geno, models = "gblup", k = 5, seed = 1))

# }
```
