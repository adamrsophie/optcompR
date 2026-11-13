# Add parameters to a model

Parameters are fixed numbers that enter expressions symbolically, so the
model can be solved again for new data with
[`opt_set_parameter()`](https://adamrsophie.github.io/optcompR/reference/opt_set_parameter.md)
without rebuilding its constraints. Problem classification and convexity
checks happen after the current values are substituted, so a parameter
that multiplies a quadratic term still gives a quadratic program.

## Usage

``` r
opt_parameter(model, name, value, dim = NULL, dimnames = NULL)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

- name:

  Block name. Must be unique among the model's variables and parameters.
  Element names are `name[i]` or `name[i,j]`, or use `dimnames`.

- value:

  Numeric values, recycled from length one.

- dim:

  Dimensions; by default the length (or
  [`dim()`](https://rdrr.io/r/base/dim.html)) of `value`.

- dimnames:

  Optional character labels, a vector for one dimension or a list with
  one entry per dimension.

## Value

An `opt_expr` referring to the parameters.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x")
p <- opt_parameter(m, "p", value = 1)
opt_objective(m, (x - p)^2)
opt_value(opt_solve(m), x)
#> [1] 1
opt_set_parameter(m, p, 5)
opt_value(opt_solve(m), x)
#> [1] 5
```
