# Set the objective of a model

Set the objective of a model

## Usage

``` r
opt_objective(model, expr, sense = NULL)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

- expr:

  A scalar expression. Constant terms are kept and reported in objective
  values.

- sense:

  Optional new sense, `"min"` or `"max"`; by default the model's sense
  is kept.

## Value

The model, invisibly.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x", lower = 1)
opt_objective(m, 2 * x + 5)
opt_objective_value(opt_solve(m))
#> [1] 7
```
