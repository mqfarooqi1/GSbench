# Genomic relationship matrix (VanRaden)

Builds the additive genomic relationship matrix from allele dosages
using VanRaden's first method: with markers centred by twice the allele
frequency, `G = W W' / (2 * sum(p (1 - p)))`. This is the standard
additive GRM used in GBLUP.

## Usage

``` r
Gmatrix(geno, min_maf = 0)
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

An `n` x `n` symmetric genomic relationship matrix with the row/column
names taken from `rownames(geno)`.

## References

VanRaden, P. M. (2008) "Efficient methods to compute genomic
predictions." Journal of Dairy Science 91, 4414-4423.
[doi:10.3168/jds.2007-0980](https://doi.org/10.3168/jds.2007-0980)

## Examples

``` r
sim <- simulate_population(n = 80, m = 400, seed = 1)
G <- Gmatrix(sim$geno)
dim(G)
#> [1] 80 80
```
