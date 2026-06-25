# Stacked super-learner ensemble of genomic prediction models

Combines several base models into one predictor by *stacking*: each base
model's out-of-fold cross-validated predictions are used to learn a set
of non-negative weights (constrained to sum to one), and the final
prediction is that weighted average of the base models refit on all the
data. This is the Breiman / van der Laan stacked-regression
(super-learner) idea applied to genomic selection; in practice it tends
to match or beat the best single model without having to know in advance
which that is.

## Usage

``` r
gs_ensemble(y, geno, base_models = NULL, inner_k = 5, seed = NULL, ...)
```

## Arguments

- y:

  Numeric phenotype vector (length n), no missing values.

- geno:

  Marker matrix (n x m, 0/1/2, no missing values).

- base_models:

  Character vector of base model names. Defaults to every available
  model except the ensemble itself.

- inner_k:

  Folds for the inner stacking cross-validation. Default 5.

- seed:

  Optional seed for the inner folds (via
  [`withr::with_seed()`](https://withr.r-lib.org/reference/with_seed.html)).

- ...:

  Passed to
  [`gs_fit()`](https://mqfarooqi1.github.io/GSbench/reference/gs_fit.md)
  for the base models.

## Value

An object of class `gs_ensemble` (and `gs_model`): a list with
`base_names`, `weights` (named, summing to 1), the refit `base_fits`,
and the out-of-fold prediction matrix `oof`.

## References

van der Laan, M. J., Polley, E. C. and Hubbard, A. E. (2007) "Super
Learner." Statistical Applications in Genetics and Molecular Biology 6,
Article 25.
[doi:10.2202/1544-6115.1309](https://doi.org/10.2202/1544-6115.1309)

## Examples

``` r
sim <- simulate_population(n = 100, m = 300, seed = 1)
ens <- gs_ensemble(sim$pheno, sim$geno, base_models = "gblup", seed = 1)
ens$weights
#> gblup 
#>     1 
```
