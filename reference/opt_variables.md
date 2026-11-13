# List the variables of a model

List the variables of a model

## Usage

``` r
opt_variables(model)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

## Value

A data frame with one row per scalar variable: name, block, type, lower
and upper bound, and start value.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x", dim = 2, lower = 0)
opt_variables(m)
#>   name block       type lower upper start
#> 1 x[1]     x continuous     0   Inf    NA
#> 2 x[2]     x continuous     0   Inf    NA
```
