# Convert an optimization expression to a symcompR expression

Builds the symbolic form of a scalar expression, with one symcompR
`Symbol` per variable or parameter (named as in the model). Coefficients
are double-precision floats unless they are integers. The symbolic form
can be simplified or differentiated with symcompR and converted back
with
[`opt_from_sym()`](https://adamrsophie.github.io/optcompR/reference/opt_from_sym.md).

## Usage

``` r
opt_to_sym(x)
```

## Arguments

- x:

  A scalar `opt_expr`.

## Value

A symcompR expression.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x")
s <- opt_to_sym(3 * x^2 + exp(x))
diff(s, symcompR::Symbol("x"))
#> 6*x + exp(x)
```
