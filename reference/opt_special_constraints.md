# List the special constraints of a model

Cones, indicator constraints and special ordered sets are stored apart
from the constraint rows listed by
[`opt_constraints()`](https://adamrsophie.github.io/optcompR/reference/opt_constraints.md).
Delete one with
[`opt_delete_constraint()`](https://adamrsophie.github.io/optcompR/reference/opt_delete_constraint.md),
passing the reference returned when it was added or one rebuilt with
`opt_special_constraint()`.

## Usage

``` r
opt_special_constraints(model)

opt_special_constraint(model, kind, name)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

- kind:

  Kind of the special constraint.

- name:

  Its name, as listed by `opt_special_constraints()`.

## Value

A data frame with one row per special constraint: `kind` (`"cone"`,
`"indicator"` or `"sos"`), `type` (for example `"soc"`, `"exp"` or
`"sos2"`), `name` and `size` (the number of scalar arguments).

For `opt_special_constraint()`, a reference usable with
[`opt_delete_constraint()`](https://adamrsophie.github.io/optcompR/reference/opt_delete_constraint.md).

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x", dim = 3, lower = 0, upper = 1)
s <- opt_sos1(m, x)
opt_special_constraints(m)
#>   kind type name size
#> 1  sos sos1 sos1    3
opt_delete_constraint(m, s)
nrow(opt_special_constraints(m))
#> [1] 0
```
