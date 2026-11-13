# Test whether an object is an optimization expression

Test whether an object is an optimization expression

## Usage

``` r
is_opt_expr(x)
```

## Arguments

- x:

  Any object.

## Value

`TRUE` if `x` is an `opt_expr`.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x")
is_opt_expr(2 * x + 1)
#> [1] TRUE
is_opt_expr(3)
#> [1] FALSE
```
