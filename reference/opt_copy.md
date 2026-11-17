# Copy a model

Models have reference semantics, so assigning a model to a new name does
not copy it. `opt_copy()` returns an independent model with the same
variables, parameters, constraints, special constraints, objective and
values; changing one model leaves the other unchanged. Expressions built
from the original still refer to the original, so retrieve the blocks of
the copy with
[`opt_get()`](https://adamrsophie.github.io/optcompR/reference/opt_get.md).

## Usage

``` r
opt_copy(model, name = NULL)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

- name:

  Optional name of the copy; by default the original name.

## Value

A new `opt_model`.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x", lower = 0)
opt_objective(m, x)
m2 <- opt_copy(m)
opt_set_bounds(m2, opt_get(m2, "x"), lower = 2)
opt_objective_value(opt_solve(m))   # 0
#> [1] 0
opt_objective_value(opt_solve(m2))  # 2
#> [1] 2
```
