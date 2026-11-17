# Inspect the terms of an expression

`opt_terms()` returns the polynomial structure of an expression: its
constant, linear and quadratic coefficients per element, using the names
of the model's variables and parameters. `opt_degree()` returns the
polynomial degree of each element (0, 1 or 2, and `NA` for elements with
nonlinear terms; parameters count like variables), and `opt_symbols()`
the names of the variables and parameters an expression depends on.

## Usage

``` r
opt_terms(x)

opt_degree(x)

opt_symbols(x)
```

## Arguments

- x:

  An `opt_expr` or a number.

## Value

For `opt_terms()`, a list with `constant` (one value per element),
`linear` (a data frame with columns `element`, `symbol` and `coef`),
`quadratic` (columns `element`, `symbol1`, `symbol2`, `coef`, meaning
`coef * symbol1 * symbol2`) and `nonlinear` (`TRUE` when the expression
has nonlinear nodes, which `opt_terms()` cannot list). For
`opt_degree()` an integer vector; for `opt_symbols()` a character
vector.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x", dim = 2)
e <- c(3 * x[1] + 2, x[1] * x[2] - x[2])
opt_terms(e)
#> $constant
#> [1] 2 0
#> 
#> $linear
#>   element symbol coef
#> 1       1   x[1]    3
#> 2       2   x[2]   -1
#> 
#> $quadratic
#>   element symbol1 symbol2 coef
#> 1       2    x[1]    x[2]    1
#> 
#> $nonlinear
#> [1] FALSE
#> 
opt_degree(c(e, exp(x[1])))
#> [1]  1  2 NA
opt_symbols(e[2])
#> [1] "x[1]" "x[2]"
```
