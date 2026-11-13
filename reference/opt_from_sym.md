# Convert a symcompR expression to an optimization expression

Translates a symbolic expression from the symcompR package into an
`opt_expr` of a model. Symbols are matched by name to the model's
variables and parameters (scalar elements such as `x[2]` included), or
to the expressions in `map`. Supported nodes are sums, products, powers,
integers, rationals, floats, numeric constants such as pi, and the
functions `exp`, `log`, `sqrt`, `sin`, `cos`, `tan`, `tanh`, `sinh`,
`cosh`, `asin`, `acos`, `atan` and `Abs`. Any other node is an error, so
a converted expression always means the same function as the symbolic
one. Rational numbers become double-precision values.

## Usage

``` r
opt_from_sym(model, expr, map = list())
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

- expr:

  A symcompR expression (class `sym`) or number.

- map:

  Optional named list from symbol names to expressions or numbers; it
  takes precedence over model names.

## Value

An `opt_expr`.

## See also

[`opt_to_sym()`](https://adamrsophie.github.io/optcompR/reference/opt_to_sym.md)

## Examples

``` r
library(symcompR)
m <- opt_model()
x <- opt_variable(m, "x")
y <- opt_variable(m, "y")
s <- (Symbol("x") - 1)^2 + Symbol("x") * Symbol("y") + exp(Symbol("y"))
e <- opt_from_sym(m, s)
opt_objective(m, e)
opt_classify(m)$class
#> [1] "NLP"
```
