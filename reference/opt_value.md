# Values of expressions at a solution

Evaluates any expression of the model (variables, slices, or derived
expressions) at the returned primal point, using the parameter values
that were in force when the model was solved. Later parameter changes do
not affect stored results.

## Usage

``` r
opt_value(result, x)
```

## Arguments

- result:

  An `opt_result`.

- x:

  An expression of the solved model, or a variable block name.

## Value

A numeric vector, or an array for matrix expressions.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x", dim = 2, lower = 1)
opt_objective(m, sum(x))
r <- opt_solve(m)
opt_value(r, x)
#> [1] 1 1
opt_value(r, 2 * x[1] + 1)
#> [1] 3
```
