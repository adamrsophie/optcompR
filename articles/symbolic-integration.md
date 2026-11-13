# Symbolic expressions with symcompR

``` r

library(optcompR)
```

optcompR describes an optimization problem once, as variables,
expressions, constraints and an objective, and then classifies it,
chooses a capable solver

symcompR is a computer algebra system for R. optcompR converts
expressions in both directions, so a model can be derived symbolically
and solved numerically, or a model expression can be simplified and
differentiated symbolically.

## From symcompR to optcompR

Symbols map to variables with the same name (or through `map`).

``` r

library(symcompR)
m <- opt_model()
x <- opt_variable(m, "x")
y <- opt_variable(m, "y")
s <- (Symbol("x") - 1)^2 + Symbol("x") * Symbol("y") + exp(Symbol("y"))
e <- opt_from_sym(m, s)
e
#> <opt_expr 1, nonlinear>
#>   ((exp(y) + (x*y + 0*x + 0*y)) + (x^2 - x - x + 1))
opt_objective(m, e)
opt_classify(m)$class
#> [1] "NLP"
```

## From optcompR to symcompR

``` r

m2 <- opt_model()
z <- opt_variable(m2, "z")
f <- 3 * z^2 + exp(z)
sf <- opt_to_sym(f)
sf
#> 3*z^2 + exp(z)
diff(sf, Symbol("z"))
#> 6*z + exp(z)
```

## Checking derivatives symbolically

The symbolic derivative and optcompR’s exact gradient agree:

``` r

opt_objective(m2, f)
opt_gradient(m2, c(z = 0.7))
#>        z 
#> 6.213753
6 * 0.7 + exp(0.7)
#> [1] 6.213753
```

External functions created with
[`opt_function()`](https://adamrsophie.github.io/optcompR/reference/opt_function.md)
have no symbolic form and cannot be converted.
