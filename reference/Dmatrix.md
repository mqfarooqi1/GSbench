# Dominance genomic relationship matrix

Builds the dominance genomic relationship matrix from allele dosages
using the genotypic parameterisation of Vitezica et al. (2013):
homozygotes and heterozygotes are coded so the matrix captures dominance
(heterozygote versus homozygote) deviations, complementing the additive
relationships from
[`Gmatrix()`](https://mqfarooqi1.github.io/GSbench/reference/Gmatrix.md).
Supplying it to
[`gblup()`](https://mqfarooqi1.github.io/GSbench/reference/gblup.md) as
`K` fits a dominance GBLUP model, and the additive and dominance GEBVs
can be compared.

## Usage

``` r
Dmatrix(geno, min_maf = 0)
```

## Arguments

- geno:

  A numeric marker matrix (individuals x markers), coded 0/1/2, with no
  missing values (run
  [`qc_markers()`](https://mqfarooqi1.github.io/GSbench/reference/qc_markers.md)
  or
  [`impute_markers()`](https://mqfarooqi1.github.io/GSbench/reference/impute_markers.md)
  first).

- min_maf:

  Markers with minor allele frequency below this are dropped before
  building the matrix. Default 0 (keep all supplied markers).

## Value

An `n` x `n` symmetric dominance relationship matrix, with row/column
names taken from `rownames(geno)`.

## References

Vitezica, Z. G., Varona, L. & Legarra, A. (2013) "On the additive and
dominant variance and covariance of individuals within the genomic
selection scope." Genetics 195, 1223-1230.
[doi:10.1534/genetics.113.155176](https://doi.org/10.1534/genetics.113.155176)

## See also

[`Gmatrix()`](https://mqfarooqi1.github.io/GSbench/reference/Gmatrix.md),
[`gblup()`](https://mqfarooqi1.github.io/GSbench/reference/gblup.md)

## Examples

``` r
sim <- simulate_population(n = 80, m = 400, seed = 1)
D <- Dmatrix(sim$geno)
dim(D)
#> [1] 80 80
# dominance GBLUP:
fit <- gblup(sim$pheno, K = D)
```
