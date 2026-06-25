# Impute missing marker genotypes

Replaces missing values in each marker (column) with that marker's mean
dosage (i.e. twice the estimated allele frequency). Simple and fast; for
production work a model-based imputation upstream is preferable.

## Usage

``` r
impute_markers(geno)
```

## Arguments

- geno:

  A numeric marker matrix (individuals x markers), 0/1/2, possibly with
  `NA`s.

## Value

The matrix with `NA`s filled by column means. Columns that are entirely
missing are filled with 0.

## Examples

``` r
g <- matrix(c(0, 1, NA, 2, 2, 0), nrow = 3)
impute_markers(g)
#>      [,1] [,2]
#> [1,]  0.0    2
#> [2,]  1.0    2
#> [3,]  0.5    0
```
