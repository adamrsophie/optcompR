# Solver comparison

``` r

library(optcompR)
```

optcompR describes an optimization problem once, as variables,
expressions, constraints and an objective, and then classifies it,
chooses a capable solver Every backend is registered with the problem
classes it accepts, the classes it solves to global optimality and the
features it supports.

``` r

opt_solvers()
#>   solver   package available                     classes duals
#> 1  highs     highs      TRUE                  LP MILP QP  TRUE
#> 2   scip      scip      TRUE LP MILP QP MIQP QCQP MIQCQP FALSE
#> 3   osqp      osqp      TRUE                       LP QP  TRUE
#> 4 nloptr    nloptr      TRUE              NLP QCQP QP LP FALSE
#> 5    roi ROI, slam      TRUE                  LP MILP QP FALSE
#>                                                                            description
#> 1       HiGHS: dual simplex, interior point and branch and bound (LP, convex QP, MILP)
#> 2 SCIP: branch and bound with spatial branching, global for nonconvex quadratic models
#> 3       OSQP: operator splitting for convex QP and LP (first-order, moderate accuracy)
#> 4                    NLopt SLSQP (default): local gradient-based nonlinear programming
#> 5       ROI: R Optimization Infrastructure with installed plugins (options$roi_solver)
opt_solver_capabilities("highs")
#> $name
#> [1] "highs"
#> 
#> $package
#> [1] "highs"
#> 
#> $available
#> [1] TRUE
#> 
#> $version
#> [1] "1.15.1.2"
#> 
#> $classes
#> [1] "LP"   "MILP" "QP"  
#> 
#> $convex_only
#> [1] "QP"
#> 
#> $global
#> [1] "LP"   "MILP" "QP"  
#> 
#> $features
#> $features$duals
#> [1] TRUE
#> 
#> $features$ranging
#> [1] TRUE
#> 
#> $features$warm_start
#> [1] TRUE
#> 
#> $features$time_limit
#> [1] TRUE
#> 
#> $features$mip_gap
#> [1] TRUE
#> 
#> $features$iteration_limit
#> [1] TRUE
#> 
#> 
#> $description
#> [1] "HiGHS: dual simplex, interior point and branch and bound (LP, convex QP, MILP)"
```

## Support matrix

| Backend | Classes | Global for | Duals | Ranging | Warm start |
|----|----|----|----|----|----|
| highs | LP, MILP, convex QP | all of these | LP, QP | LP | yes |
| scip | LP, MILP, QP, MIQP, QCQP, MIQCQP | all, including nonconvex | no | no | no |
| osqp | convex QP, LP | convex problems | yes | no | no |
| nloptr | NLP, QCQP, QP, LP (SLSQP) | convex problems | estimated by [`opt_check()`](https://adamrsophie.github.io/optcompR/reference/opt_check.md) | no | start point |
| roi | whatever the chosen ROI plugin solves | plugin | plugin | no | no |

Not supported in this version: conic programs (SOCP, SDP), general
MINLP, indicator and SOS constraints, IPOPT.

## Automatic choice

`solver = "auto"` picks the highest priority backend that accepts the
class and convexity of the model. For nonconvex problems global solvers
come first. When nothing fits, the error explains why each backend was
rejected:

``` r

m <- opt_model()
x <- opt_variable(m, "x", type = "integer", lower = 0, upper = 3)
opt_objective(m, exp(x))
opt_solve(m)
#> Error:
#> ! No installed solver supports this MINLP (unknown). highs does not support MINLP problems; scip does not support MINLP problems; osqp does not support MINLP problems; nloptr does not support MINLP problems; roi does not support MINLP problems.
```

## Comparing backends on one model

``` r

m <- opt_model()
x <- opt_variable(m, "x", dim = 3, lower = 0)
opt_constraint(m, sum(x) >= 1)
opt_objective(m, sum((x - c(0.2, 0.5, 0.1))^2) + x[1] * x[2])
solvers <- intersect(c("highs", "osqp", "scip", "nloptr"), with(opt_solvers(), solver[available]))
res <- lapply(solvers, function(s) opt_solve(m, solver = s))
data.frame(solver = solvers,
           status = vapply(res, opt_status, ""),
           objective = vapply(res, opt_objective_value, 0),
           max_row_violation = vapply(res, function(r) opt_check(r)$max_row_violation, 0))
#>   solver  status objective max_row_violation
#> 1  highs optimal 0.1171429      1.110223e-16
#> 2   osqp optimal 0.1171429      0.000000e+00
#> 3   scip optimal 0.1171429      0.000000e+00
#> 4 nloptr optimal 0.1171429      0.000000e+00
```

The objective in every result is evaluated from the model, not copied
from the solver, so the numbers are comparable.
`opt_objective_value(r, solver = TRUE)` returns the solver’s own value.

## ROI plugins

``` r

m <- opt_model("max")
x <- opt_variable(m, "x", dim = 2, type = "integer", lower = 0)
opt_constraint(m, 2 * x[1] + 3 * x[2] <= 12)
opt_objective(m, 3 * x[1] + 4 * x[2])
r <- opt_solve(m, solver = "roi", options = list(roi_solver = "glpk"))
opt_value(r, x)
#> [1] 6 0
```

## Solver options

Options are passed through to the backend; for HiGHS they override
defaults:

``` r

r <- opt_solve(m, solver = "highs", options = list(presolve = "off"))
opt_status(r)
#> [1] "optimal"
```
