# Add constraints to a model

Adds a block of constraints. The block is either a relation built with
`<=`, `>=` or `==` (for example `A %*% x <= b`), or an expression with
explicit `lower` and `upper` bounds, which gives ranged rows
`lower <= expr <= upper`. Vector relations add one row per element.

## Usage

``` r
opt_constraint(model, relation, name = NULL, lower = NULL, upper = NULL)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

- relation:

  An `opt_relation` from a comparison, or an `opt_expr` when `lower` or
  `upper` is given.

- name:

  Block name, unique within the model. Defaults to `c1`, `c2`, and so
  on.

- lower, upper:

  Row bounds for the ranged form, recycled from length one.

## Value

An `opt_constraint` reference, used by
[`opt_dual()`](https://adamrsophie.github.io/optcompR/reference/opt_dual.md),
[`opt_slack()`](https://adamrsophie.github.io/optcompR/reference/opt_slack.md)
and
[`opt_delete_constraint()`](https://adamrsophie.github.io/optcompR/reference/opt_delete_constraint.md).
Returned invisibly.

## Details

Constraints are stored exactly as written. They are never approximated;
a nonlinear constraint makes the model nonlinear, and solvers that
cannot handle it are not selected.

## Examples

``` r
m <- opt_model("max")
x <- opt_variable(m, "x", dim = 2, lower = 0)
cap <- opt_constraint(m, x[1] + x[2] <= 4, name = "capacity")
opt_constraint(m, x[1] - x[2], lower = -1, upper = 1, name = "balance")
opt_objective(m, 3 * x[1] + 2 * x[2])
r <- opt_solve(m)
opt_dual(r, cap)
#> capacity 
#>      2.5 
```
