# Solver status of a result

`opt_status()` returns the termination status, one of:

## Usage

``` r
opt_status(result)

opt_primal_status(result)

opt_dual_status(result)

opt_has_solution(result)
```

## Arguments

- result:

  An `opt_result`.

## Value

A string.

## Details

- `"optimal"`: the solver proved optimality within its tolerances
  (globally, for the classes the backend solves globally);

- `"locally_optimal"`: a local solver converged; nonconvex problems may
  have better solutions elsewhere;

- `"almost_optimal"`: the solver met relaxed tolerances only (OSQP's
  "solved inaccurate");

- `"locally_infeasible"`: a local solver stopped at a point that
  violates constraints; this is not a proof of infeasibility;

- `"infeasible"`, `"unbounded"`, `"infeasible_or_unbounded"`: as
  certified by the solver; HiGHS may report infeasible or unbounded
  without deciding which;

- `"time_limit"`, `"iteration_limit"`, `"node_limit"`,
  `"solution_limit"`, `"objective_limit"`, `"memory_limit"`,
  `"interrupted"`: the solver stopped early; a feasible point may still
  exist, see `opt_primal_status()`;

- `"numerical_error"`, `"invalid_model"`, `"unknown"`.

A mixed-integer incumbent found before a limit is never reported as
optimal. `opt_primal_status()` and `opt_dual_status()` describe the
returned points: `"feasible_point"`, `"infeasible_point"`,
`"no_solution"`, `"not_available"` or `"unknown"`.

## Examples

``` r
m <- opt_model()
x <- opt_variable(m, "x", lower = 0)
opt_constraint(m, x <= -1)
opt_objective(m, x)
opt_status(opt_solve(m))
#> [1] "infeasible"
```
