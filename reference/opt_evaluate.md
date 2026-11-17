# Evaluate an expression at given values

Evaluate an expression at given values

## Usage

``` r
opt_evaluate(x, values)
```

## Arguments

- x:

  An `opt_expr`.

- values:

  Values of the model's variables: a numeric vector with one value per
  variable (in the order of
  [`opt_variables()`](https://adamrsophie.github.io/optcompR/reference/opt_variables.md))
  or a named vector naming every variable. Parameters take their current
  values.

## Value

A numeric vector (or matrix, for matrix expressions).

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x", dim = 2)
opt_evaluate(x[1]^2 + 3 * x[2], c(2, 1))
#> [1] 7
```
