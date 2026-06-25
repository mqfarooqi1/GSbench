# Fit GBLUP by REML

Fits the genomic best linear unbiased prediction model
`y = X b + g + e`, with `g ~ N(0, G s2u)` and `e ~ N(0, I s2e)`,
estimating the variance components by restricted maximum likelihood
(REML). The solver uses the spectral / eigendecomposition approach of
Endelman (2011) (the same method as
[`rrBLUP::mixed.solve`](https://rdrr.io/pkg/rrBLUP/man/mixed.solve.html)):
the REML log-likelihood is profiled to a one-dimensional optimisation
over the variance ratio `lambda = s2e / s2u`.

## Usage

``` r
gblup(y, geno = NULL, K = NULL, X = NULL)
```

## Arguments

- y:

  Numeric phenotype vector (length n), no missing values.

- geno:

  Optional marker matrix (n x m, coded 0/1/2). Used to build `K` and to
  derive marker effects for out-of-sample prediction.

- K:

  Optional n x n genomic relationship matrix. If `NULL`, built from
  `geno`. If you pass `K` directly (no `geno`), out-of-sample prediction
  is not available.

- X:

  Optional fixed-effects design matrix (n x p). Defaults to an
  intercept.

## Value

An object of class `gblup`: a list with `Vu`, `Ve`, `h2`, `beta` (fixed
effects), `gebv` (length-n GEBVs), `lambda`, `K`, and — when `geno` was
given — `marker_effects`, `marker_means` (2p) and `intercept`.

## Details

When `geno` is supplied, the equivalent ridge-regression marker effects
are also returned, so the fit can predict breeding values for *new*
genotypes via
[`predict.gblup()`](https://mqfarooqi1.github.io/GSbench/reference/predict.gblup.md).

## References

Endelman, J. B. (2011) "Ridge regression and other kernels for genomic
selection with R package rrBLUP." The Plant Genome 4, 250-255.
[doi:10.3835/plantgenome2011.08.0024](https://doi.org/10.3835/plantgenome2011.08.0024)

## Examples

``` r
sim <- simulate_population(n = 120, m = 500, seed = 1)
fit <- gblup(sim$pheno, sim$geno)
fit$h2
#> [1] 0.6814213
cor(fit$gebv, sim$bv)
#> [1] 0.773269
```
