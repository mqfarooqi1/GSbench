# GSbench: benchmarking genomic selection and machine-learning models

GSbench provides a single, consistent interface for genomic prediction
from SNP markers. It implements GBLUP and ridge-regression BLUP in base
R, wraps machine-learning predictors (elastic net, random forest,
gradient boosting) behind the same interface, and adds a stacked
ensemble. Models are compared with breeding-relevant cross-validation
and honest accuracy reporting.

## Where to start

[`simulate_population()`](https://mqfarooqi1.github.io/GSbench/reference/simulate_population.md)
makes a toy dataset,
[`gblup()`](https://mqfarooqi1.github.io/GSbench/reference/gblup.md)
fits the GBLUP baseline, and (in later sections of the package) the
cross-validation and benchmarking functions compare models.

## Design

The package uses the S3 object system. Core numerical methods (the
genomic relationship matrix and the mixed-model solver) are written in
base R with no compiled code, so the package installs without a
toolchain. Machine-learning backends are optional (declared under
Suggests) and used only if installed.

## See also

Useful links:

- <https://github.com/mqfarooqi1/GSbench>

- Report bugs at <https://github.com/mqfarooqi1/GSbench/issues>

## Author

**Maintainer**: Muhammad Farooqi <mqfarooqi@gmail.com>

Authors:

- Muhammad Farooqi <mqfarooqi@gmail.com>
