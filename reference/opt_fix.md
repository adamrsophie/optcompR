# Fix variables to values

Sets both bounds to `value`.

## Usage

``` r
opt_fix(model, x, value)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

- x:

  A variable expression (or slice) or the variable block name.

- value:

  Values, recycled from length one.

## Value

The model, invisibly.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x", lower = 0)
opt_fix(m, x, 3)
opt_variables(m)
#>   name block       type lower upper start
#> 1    x     x continuous     3     3    NA
```
