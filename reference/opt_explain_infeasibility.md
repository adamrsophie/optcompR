# Explain why a model is infeasible

Solves an elastic relaxation of the model: every constraint row gets
nonnegative slack variables that let it move toward feasibility, and
with `bounds = TRUE` the bounds of continuous and integer variables do
too. The relaxation minimizes the weighted total (L1) violation, so its
solution names rows and bounds that can move, and by how much, for the
model to become feasible at the least weighted total change. This is not
the smallest number of rows to change, nor an irreducible infeasible
subset: two copies of `x >= 1` against one `100 * x <= 0` relax the two
copies (total change 2) although relaxing the single upper row would
suffice (change 100). Special constraints (indicators and SOS), binary
bounds and integrality are kept as they are.

## Usage

``` r
opt_explain_infeasibility(
  model,
  bounds = TRUE,
  weights = NULL,
  solver = "auto",
  ...,
  tol = 1e-07
)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

- bounds:

  Also relax variable bounds.

- weights:

  Optional named numeric vector of weights by constraint block name (and
  `".bounds"` for variable bounds); the default weight is

  1.  A larger weight makes moving that block less attractive.

- solver, ...:

  Passed to
  [`opt_solve()`](https://adamrsophie.github.io/optcompR/reference/opt_solve.md)
  for the relaxation.

- tol:

  Violations at or below `tol` are reported as zero.

## Value

An `opt_infeasibility` list with `feasible` (whether the model is
feasible, that is, no unweighted violation above `tol` is needed; the
weights do not change this verdict), `total` (the weighted sum of
violations), `violation` (their unweighted sum), `violations`, a data
frame with columns `kind` (`"row"` or `"bound"`), `name`, `block`,
`side` (`"lower"` or `"upper"`), `bound` (the original bound), `value`
(the row or variable value at the relaxed solution) and `violation`, and
`result`, the
[`opt_solve()`](https://adamrsophie.github.io/optcompR/reference/opt_solve.md)
result of the relaxation.

## Details

The model itself is not changed; the relaxation is built on a copy (see
[`opt_copy()`](https://adamrsophie.github.io/optcompR/reference/opt_copy.md)).
For nonlinear models the relaxation is solved by a local solver unless
the model is proven convex, so the explanation is then a local one.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x", dim = 2, lower = 0)
opt_constraint(m, x[1] + x[2] <= 1, name = "capacity")
opt_constraint(m, x[1] >= 2, name = "demand_1")
opt_constraint(m, x[2] >= 0.5, name = "demand_2")
opt_explain_infeasibility(m)
#> <opt_infeasibility: 1 row(s) and 0 bound(s) must move, total violation 1.5>
#>   row demand_1: lower bound 2, relaxed value 0.5 (violation 1.5)
```
