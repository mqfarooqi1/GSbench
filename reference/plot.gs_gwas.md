# Plot a GWAS scan (Manhattan or QQ)

Plot a GWAS scan (Manhattan or QQ)

## Usage

``` r
# S3 method for class 'gs_gwas'
plot(x, type = c("manhattan", "qq"), alpha = 0.05, ...)
```

## Arguments

- x:

  A `gs_gwas` object from
  [`gwas()`](https://mqfarooqi1.github.io/GSbench/reference/gwas.md).

- type:

  `"manhattan"` (default) or `"qq"`.

- alpha:

  Genome-wide significance level for the Bonferroni line on the
  Manhattan plot (default 0.05).

- ...:

  Further graphical parameters passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

The `gs_gwas` object, invisibly.

## Examples

``` r
sim <- simulate_population(n = 150, m = 500, n_qtl = 10, h2 = 0.6, seed = 1)
scan <- gwas(sim$pheno, sim$geno, n_pc = 3)
plot(scan)

plot(scan, type = "qq")
```
