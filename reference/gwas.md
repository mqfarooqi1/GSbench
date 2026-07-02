# Single-marker genome-wide association scan

Tests each marker for additive association with a phenotype, optionally
correcting for population structure with genotype principal components.
This complements the whole-genome prediction models in GSbench:
prediction asks *how well can I predict*, whereas GWAS asks *which
markers are associated*.

## Usage

``` r
gwas(y, geno, n_pc = 0, covariates = NULL, map = NULL, min_maf = 0)

# S3 method for class 'gs_gwas'
print(x, top = 6, ...)
```

## Arguments

- y:

  Numeric phenotype vector (length n). Missing values are dropped.

- geno:

  Marker matrix (n x m), coded 0/1/2, with no missing values.

- n_pc:

  Number of genotype principal components to include as covariates for
  structure correction (default 0).

- covariates:

  Optional numeric matrix of additional fixed covariates (n x k).

- map:

  Optional data frame of per-marker annotation; columns `chrom` and
  `pos` (if present) are carried through and enable chromosome-aware
  Manhattan plots.

- min_maf:

  Drop markers with minor-allele frequency below this before testing
  (default 0).

- x:

  A `gs_gwas` object.

- top:

  Number of top associations to show (default 6).

- ...:

  Further arguments passed to plotting functions.

## Value

An object of class `gs_gwas` (a data frame) with one row per marker:
`marker`, optionally `chrom` and `pos`, `effect`, `se`, `statistic`,
`p_value` and `log10p`. It has
[`print()`](https://rdrr.io/r/base/print.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods
(Manhattan and QQ plots); see
[`plot.gs_gwas()`](https://mqfarooqi1.github.io/GSbench/reference/plot.gs_gwas.md).

## Details

For speed the phenotype and the markers are projected onto the residual
space of the covariates (intercept, optional principal components and
extra covariates), and the partial regression slope of each marker is
tested with a t-test. This is algebraically equivalent to fitting one
linear model per marker but avoids the per-marker loop.

## Functions

- `print(gs_gwas)`: Print a short summary and the top associations.

## References

Yu, J. et al. (2006) "A unified mixed-model method for association
mapping that accounts for multiple levels of relatedness." Nature
Genetics 38, 203-208.
[doi:10.1038/ng1702](https://doi.org/10.1038/ng1702)

## See also

[`gblup()`](https://mqfarooqi1.github.io/GSbench/reference/gblup.md) for
whole-genome prediction;
[`plot.gs_gwas()`](https://mqfarooqi1.github.io/GSbench/reference/plot.gs_gwas.md).

## Examples

``` r
sim <- simulate_population(n = 150, m = 500, n_qtl = 10, h2 = 0.6, seed = 1)
scan <- gwas(sim$pheno, sim$geno, n_pc = 3)
head(scan[order(scan$p_value), ])
#> <gs_gwas: 6 markers, n = 150>
#>   6 markers pass a Bonferroni threshold (p < 0.0083)
#>   top associations:
#>  marker    effect        se statistic      p_value   log10p
#>  snp138 -2.308998 0.4529788 -5.097364 1.058122e-06 5.975464
#>  snp131 -2.013127 0.4696586 -4.286363 3.294649e-05 4.482191
#>  snp470 -2.083612 0.4866277 -4.281738 3.356011e-05 4.474177
#>   snp41 -1.996625 0.5136541 -3.887101 1.538910e-04 3.812787
#>  snp222 -2.142854 0.5927173 -3.615305 4.131551e-04 3.383887
#>  snp480  1.530068 0.4465508  3.426413 7.959736e-04 3.099101
```
