# Simulate a genomic prediction dataset

Generates a small marker matrix and a phenotype with a known
narrow-sense heritability, for examples, tests and demonstrations.
Markers are coded as allele dosages (0, 1, 2). A random subset of
markers are QTL with additive effects; the phenotype is the resulting
breeding value plus normal noise scaled to the target heritability.

## Usage

``` r
simulate_population(n = 200, m = 1000, n_qtl = 50, h2 = 0.5, seed = NULL)
```

## Arguments

- n:

  Number of individuals. Default 200.

- m:

  Number of markers. Default 1000.

- n_qtl:

  Number of markers with a non-zero (QTL) effect. Default 50.

- h2:

  Target narrow-sense heritability in (0, 1\]. Default 0.5.

- seed:

  Optional integer seed for reproducibility. Applied with
  [`withr::with_seed()`](https://withr.r-lib.org/reference/with_seed.html),
  which scopes the seed locally and leaves the caller's random-number
  state unchanged.

## Value

A list with `geno` (an `n` x `m` 0/1/2 matrix, row/column names set),
`pheno` (length-`n` numeric), `bv` (true breeding values), `qtl`
(indices of the causal markers) and `h2` (the realised heritability).

## Examples

``` r
sim <- simulate_population(n = 100, m = 300, seed = 1)
dim(sim$geno)
#> [1] 100 300
length(sim$pheno)
#> [1] 100
```
