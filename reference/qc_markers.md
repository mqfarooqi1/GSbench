# Quality-control filter for marker data

Drops markers failing a call-rate or minor-allele-frequency threshold,
and (optionally) monomorphic markers, then imputes any remaining missing
values.

## Usage

``` r
qc_markers(geno, maf = 0.05, max_missing = 0.1, impute = TRUE)
```

## Arguments

- geno:

  A numeric marker matrix (individuals x markers), coded 0/1/2.

- maf:

  Minimum minor allele frequency to keep a marker. Default 0.05.

- max_missing:

  Maximum fraction of missing calls to keep a marker. Default 0.1.

- impute:

  Whether to mean-impute remaining missing values. Default `TRUE`.

## Value

A list with `geno` (the filtered, optionally imputed matrix) and
`removed` (a named integer vector counting markers dropped by each
rule).

## Examples

``` r
sim <- simulate_population(n = 50, m = 200, seed = 1)
qc <- qc_markers(sim$geno)
dim(qc$geno)
#> [1]  50 199
qc$removed
#>   call_rate         maf monomorphic 
#>           0           1           0 
```
