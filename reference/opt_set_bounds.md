# Change variable bounds

Change variable bounds

## Usage

``` r
opt_set_bounds(model, x, lower = NULL, upper = NULL)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

- x:

  A variable expression (or slice) or the variable block name.

- lower, upper:

  New bounds, recycled from length one; `NULL` keeps the current bound.

## Value

The model, invisibly.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x", dim = 2)
opt_set_bounds(m, x, lower = 0, upper = c(1, 2))
opt_variables(m)
#>   name block       type lower upper start
#> 1 x[1]     x continuous     0     1    NA
#> 2 x[2]     x continuous     0     2    NA
```
