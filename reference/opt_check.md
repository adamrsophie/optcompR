# Independent feasibility and optimality diagnostics

Recomputes, from the original model expressions and not from the solver
or the standard form, the constraint values, bound and integrality
violations and the objective value at the returned point. For continuous
problems with duals it also checks the Karush Kuhn Tucker conditions
under the sign convention documented in
[`opt_dual()`](https://adamrsophie.github.io/optcompR/reference/opt_dual.md):

## Usage

``` r
opt_check(result, tol = 1e-06)
```

## Arguments

- result:

  An `opt_result`.

- tol:

  Feasibility and optimality tolerance.

## Value

An `opt_diagnostics` list with `max_row_violation`,
`max_bound_violation`, `max_integrality_violation`, `objective_error`,
`kkt` (or `NULL`), `feasible` and `rows`, a data frame of row residuals.

## Details

- primal feasibility: rows and bounds hold within `tol`;

- stationarity: \\\nabla f(x) - J(x)^\top y - z = 0\\, with exact
  gradients;

- dual feasibility: each dual has the sign its active bound allows;

- complementary slackness: \\y_i\\ times the row slack is zero.

Violations are scaled by `max(1, |bound|)` for rows and by
`max(1, |gradient|)` for stationarity. Feasibility passes at `tol`; the
KKT test passes at `sqrt(tol)`, because duals and active sets carry the
square root of solver tolerances, and rows within `sqrt(tol)` of a bound
count as active. When a solver returns no multipliers, they are
estimated by least squares on the active set and the report says so. KKT
conditions certify global optimality only for convex problems; for
nonconvex problems they describe a stationary point.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x")
y <- opt_variable(m, "y")
opt_constraint(m, x + y >= 4)
opt_objective(m, (x - 1)^2 + (y - 2)^2)
opt_check(opt_solve(m))
#> <opt_diagnostics: feasible>
#>   max row violation (scaled):   0
#>   max bound violation (scaled): 0
#>   max integrality violation:    0
#>   objective recomputation error: 0
#>   KKT: stationarity 2.5e-07, dual feasibility 0, complementarity 0 (pass)
```
