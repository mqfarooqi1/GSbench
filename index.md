# GSbench

Benchmark genomic-selection models — classic and machine-learning — from
SNP marker data, through **one interface**, with **breeding-relevant
cross-validation** and **honest accuracy reporting**.

The problem GSbench addresses: people increasingly throw `glmnet`,
`ranger`, or `xgboost` at marker matrices, but hand-roll the
cross-validation (often incorrectly) and compare models on unequal
footing. GSbench fits the standard baselines (GBLUP, ridge marker
effects) **and** the ML methods behind a single
[`gs_fit()`](https://mqfarooqi1.github.io/GSbench/reference/gs_fit.md)/[`predict()`](https://rdrr.io/r/stats/predict.html)
API, runs them through the same CV, and reports predictive ability you
can actually trust — plus a stacked ensemble that combines them.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("mqfarooqi1/GSbench")
```

Only `graphics`, `stats` and `withr` are required. The ML backends —
`glmnet`, `ranger`, `xgboost` — are optional (Suggests); install
whichever you want to use.

## Quick start

``` r

library(GSbench)

sim <- simulate_population(n = 300, m = 2000, h2 = 0.5, seed = 1)

# one model
fit <- gs_fit(sim$pheno, sim$geno, model = "gblup")
gebv <- predict(fit, sim$geno)

# compare every available model (incl. the stacked ensemble) under one CV
bench <- gs_benchmark(sim$pheno, sim$geno, k = 5, seed = 1)
bench
plot(bench)
```

             model  mean    sd n_folds
       elastic_net 0.367 0.187       5
             gblup 0.334 0.189       5
          ensemble 0.328 0.165       5
     random_forest 0.269 0.185       5
           xgboost 0.185 0.318       5
      (accuracy = predictive ability, cor(pred, observed) on held-out data)

## What’s in it

**Core (base R, no compiled code, no heavy deps):**

| Function | Purpose |
|----|----|
| [`simulate_population()`](https://mqfarooqi1.github.io/GSbench/reference/simulate_population.md) | Reproducible SNP + phenotype simulator with known h² |
| [`qc_markers()`](https://mqfarooqi1.github.io/GSbench/reference/qc_markers.md), [`impute_markers()`](https://mqfarooqi1.github.io/GSbench/reference/impute_markers.md) | Call-rate / MAF / monomorphic filtering, mean imputation |
| [`Gmatrix()`](https://mqfarooqi1.github.io/GSbench/reference/Gmatrix.md) | VanRaden additive genomic relationship matrix |
| [`gblup()`](https://mqfarooqi1.github.io/GSbench/reference/gblup.md) | GBLUP by REML — **validated to match [`rrBLUP::mixed.solve`](https://rdrr.io/pkg/rrBLUP/man/mixed.solve.html) to 6×10⁻⁵** |

**Modelling & evaluation:**

| Function | Purpose |
|----|----|
| [`gs_fit()`](https://mqfarooqi1.github.io/GSbench/reference/gs_fit.md) / [`predict()`](https://rdrr.io/r/stats/predict.html) | Unified interface: `"gblup"`, `"elastic_net"`, `"random_forest"`, `"xgboost"`, `"ensemble"` |
| [`gs_cv()`](https://mqfarooqi1.github.io/GSbench/reference/gs_cv.md) | Cross-validation: random k-fold (CV1) or leave-one-group-out (family/environment) |
| [`gs_ensemble()`](https://mqfarooqi1.github.io/GSbench/reference/gs_ensemble.md) | **Stacked super-learner** — combines base models with non-negative CV-learned weights |
| [`gs_benchmark()`](https://mqfarooqi1.github.io/GSbench/reference/gs_benchmark.md) + [`plot()`](https://rdrr.io/r/graphics/plot.default.html) | Run all available models through one CV and compare |
| [`available_models()`](https://mqfarooqi1.github.io/GSbench/reference/available_models.md) | Which models are usable in your session |

## Why the methods are trustworthy

- **GBLUP is built from scratch in base R** (spectral REML, the Endelman
  2011 / EMMA method) and is **numerically validated against `rrBLUP`**
  in the test suite — same variance components, GEBVs correlating at
  1.0.
- **Cross-validation is the part people get wrong**, so it’s the part
  GSbench is opinionated about: correct fold construction,
  leave-group-out for family/environment structure, and accuracy
  aggregated across folds.
- The **stacked ensemble** is the Breiman / van der Laan super-learner:
  base models are combined by weights fit to their out-of-fold
  predictions (non-negative, summing to one). It tends to match or beat
  the best single model without you having to know which that is in
  advance.

## Honest limitations

- **Single trait, single environment.** Multi-trait and GxE (CV2) models
  are not here yet — that’s the obvious next direction.
- **Pure-R performance.** The GBLUP solver eigendecomposes an n×n
  matrix; fine for typical breeding populations (hundreds–few thousand
  lines), but very large panels would want a C++ backend.
- **Imputation is simple** (marker means); model-based imputation
  upstream is better for real data.
- The simulator is for demos/tests — bring your own genotypes and
  phenotypes for real work.

## References

- VanRaden, P. M. (2008) *J. Dairy Sci.* 91:4414–4423.
  <doi:10.3168/jds.2007-0980>
- Endelman, J. B. (2011) *Plant Genome* 4:250–255.
  <doi:10.3835/plantgenome2011.08.0024>
- Meuwissen, Hayes & Goddard (2001) *Genetics* 157:1819–1829.
  <doi:10.1093/genetics/157.4.1819>
- van der Laan, Polley & Hubbard (2007) *Stat. Appl. Genet. Mol. Biol.*
  6:Art.25. <doi:10.2202/1544-6115.1309>

------------------------------------------------------------------------

Muhammad Farooqi · <https://github.com/mqfarooqi1>
