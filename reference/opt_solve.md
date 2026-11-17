# Solve an optimization model

Lowers the model to its standard form, classifies it, selects a solver
that supports the class, and returns a standardized result. With
`solver = "auto"` the first installed backend (in registry priority
order) whose R binding supports the model class is used; for nonconvex
classes backends with global guarantees come first. A model is never
transformed into an approximation to fit a solver: if no installed
backend supports it, an error lists what is missing.

## Usage

``` r
opt_solve(
  model,
  solver = "auto",
  time_limit = Inf,
  mip_gap = NULL,
  verbose = FALSE,
  warm_start = TRUE,
  options = list()
)
```

## Arguments

- model:

  An
  [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md).

- solver:

  Backend name from
  [`opt_solvers()`](https://adamrsophie.github.io/optcompR/reference/opt_solvers.md),
  or `"auto"`.

- time_limit:

  Wall-clock limit in seconds, for solvers that accept it.

- mip_gap:

  Relative optimality gap for mixed-integer solvers. `NULL` keeps the
  solver default. A solve that stops with a gap above zero under a gap
  the user asked for (here or through a native option such as HiGHS
  `mip_rel_gap` or SCIP `limits/gap`), or with a gap above 1e-4 in any
  case, has status `"gap_limit"`, not `"optimal"`.

- verbose:

  Print the solver log.

- warm_start:

  Pass start values set with
  [`opt_set_start()`](https://adamrsophie.github.io/optcompR/reference/opt_set_start.md)
  or `opt_variable(start = )` to solvers that accept them, and reuse the
  native HiGHS solver and OSQP workspace of the previous solve of the
  same model: when only costs, bounds or matrix values changed (for
  example after
  [`opt_set_parameter()`](https://adamrsophie.github.io/optcompR/reference/opt_set_parameter.md)),
  they are updated in place and the solver restarts from its last basis
  or solution. `FALSE` builds a new native solver.

- options:

  Named list of backend options passed through unchanged (for HiGHS, any
  option of
  [`highs::highs_control()`](https://rdrr.io/pkg/highs/man/highs_control.html)).
  HiGHS uses one thread by default; its thread scheduler is global to
  the process, so the `threads` value of the first HiGHS solve holds for
  the whole R session and a different value is an error.

## Value

An `opt_result`. See
[`opt_status()`](https://adamrsophie.github.io/optcompR/reference/opt_status.md),
[`opt_value()`](https://adamrsophie.github.io/optcompR/reference/opt_value.md),
[`opt_objective_value()`](https://adamrsophie.github.io/optcompR/reference/opt_objective_value.md),
[`opt_dual()`](https://adamrsophie.github.io/optcompR/reference/opt_dual.md)
and
[`opt_check()`](https://adamrsophie.github.io/optcompR/reference/opt_check.md).

## Checking the returned point

The point a backend returns is checked against the model (rows from the
original expressions, bounds, integrality, cones, indicators and SOS
constraints), each side scaled by `max(1, |bound|)`, with the backend's
feasibility tolerance (1e-5 by default, 1e-4 for IPOPT and 1e-2 for the
first-order solvers OSQP and SCS). A point that fails is labelled
`"infeasible_point"`, a claimed optimum then becomes
`"locally_infeasible"` (local solvers) or `"numerical_error"`, and a
warning gives the violation. A backend is treated as local for the
classes it does not declare global unless the model is proven convex.

## Examples

``` r
m <- opt_model("max")
x <- opt_variable(m, "x", lower = 0, upper = 2)
y <- opt_variable(m, "y", lower = 0, upper = 3)
opt_constraint(m, x + y <= 4)
opt_objective(m, 3 * x + 2 * y)
r <- opt_solve(m)
r
#> <opt_result>
#>   solver:      highs 1.15.1.2
#>   problem:     LP (convex)
#>   status:      optimal (primal: feasible_point, dual: feasible_point)
#>   objective:   10
#>   time:        0.004s (build 0.002s)
opt_value(r, x)
#> [1] 2
```
