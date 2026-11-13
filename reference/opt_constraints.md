# List the constraint blocks of a model

List the constraint blocks of a model

## Usage

``` r
opt_constraints(model)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

## Value

A data frame with one row per constraint block: name, number of rows,
and the polynomial degree or nonlinear kind of its expression.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x", dim = 3, lower = 0)
opt_constraint(m, sum(x) <= 1, name = "budget")
opt_constraints(m)
#>     name rows   kind
#> 1 budget    1 affine
```
