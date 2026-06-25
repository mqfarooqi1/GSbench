# Models available in this session

Returns the genomic-prediction models GSbench can currently run.
`"gblup"` is always available; the machine-learning models require their
(suggested) package to be installed.

## Usage

``` r
available_models()
```

## Value

A character vector of usable model names.

## Examples

``` r
available_models()
#> [1] "gblup"         "elastic_net"   "random_forest" "xgboost"      
#> [5] "ensemble"     
```
