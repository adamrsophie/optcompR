# Change parameter values

Updates parameter values in place. Cached standard forms built for the
old values are invalidated; the next
[`opt_solve()`](https://adamrsophie.github.io/optcompR/reference/opt_solve.md)
rebuilds the numerical data while reusing the model structure.

## Usage

``` r
opt_set_parameter(model, x, value)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

- x:

  A parameter expression returned by
  [`opt_parameter()`](https://adamrsophie.github.io/optcompR/reference/opt_parameter.md)
  (or a slice of it), or the parameter block name.

- value:

  New values, recycled from length one.

## Value

The model, invisibly.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x", lower = 0)
b <- opt_parameter(m, "b", value = 2)
opt_constraint(m, x >= b)
opt_objective(m, x)
opt_objective_value(opt_solve(m))
#> [1] 2
opt_set_parameter(m, b, 7)
opt_objective_value(opt_solve(m))
#> [1] 7
```
