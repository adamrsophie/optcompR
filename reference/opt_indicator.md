# Indicator constraints

`opt_indicator(model, z, relation)` imposes the linear `relation` only
when the binary variable `z` equals `active`; otherwise the relation is
free. SCIP handles indicator constraints natively. For other
mixed-integer solvers they are reformulated exactly with big-M rows,
where M is computed from the variable bounds; that requires finite
bounds on the variables of the relation, and the solve stops with an
error otherwise. No approximation is made.

## Usage

``` r
opt_indicator(model, z, relation, active = 1, name = NULL)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

- z:

  A binary variable (scalar) or a vector of binaries, one per element of
  `relation`.

- relation:

  A linear relation such as `x + y <= 4`.

- active:

  `1` (the relation holds when `z = 1`) or `0`.

- name:

  Optional name.

## Value

An `opt_special_constraint` reference, invisibly.

## Examples

``` r
# a fixed cost of 10 is paid when production x is positive
m <- opt_model()
x <- opt_variable(m, "x", lower = 0, upper = 8)
on <- opt_variable(m, "on", type = "binary")
opt_indicator(m, on, x <= 0, active = 0)
opt_constraint(m, x >= 3)
opt_objective(m, 2 * x + 10 * on)
opt_value(opt_solve(m), c(x, on))
#> [1] 3 1
```
