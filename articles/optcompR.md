# Getting started with optcompR

``` r

library(optcompR)
```

optcompR describes an optimization problem once, as variables,
expressions, constraints and an objective, and then classifies it,
chooses a capable solver and returns results with a uniform meaning.
This article walks through the whole cycle on a small production
problem.

## A first model

A workshop makes tables and chairs. A table earns 7 and a chair 5.
Tables use 4 hours of carpentry and 2 of painting; chairs use 3 and 1.
There are 240 carpentry hours and 100 painting hours.

``` r

m <- opt_model("max", name = "workshop")
tables <- opt_variable(m, "tables", lower = 0)
chairs <- opt_variable(m, "chairs", lower = 0)
carpentry <- opt_constraint(m, 4 * tables + 3 * chairs <= 240, name = "carpentry")
painting <- opt_constraint(m, 2 * tables + chairs <= 100, name = "painting")
opt_objective(m, 7 * tables + 5 * chairs)
m
#> <opt_model 'workshop'>
#>   sense:       max
#>   variables:   2 (0 integer or binary)
#>   constraints: 2 rows in 2 blocks
#>   objective:   affine
```

Models have reference semantics:
[`opt_variable()`](https://adamrsophie.github.io/optcompR/reference/opt_variable.md),
[`opt_constraint()`](https://adamrsophie.github.io/optcompR/reference/opt_constraint.md)
and
[`opt_objective()`](https://adamrsophie.github.io/optcompR/reference/opt_objective.md)
change `m` in place, so there is no need to reassign it.

## Classification and solver choice

``` r

opt_classify(m)
#> <opt_classification: LP, convex>
#>   - syntactic: objective is linear; constraint rows: linear = 2
#>   - structural: linear objective
#>   - structural: all rows linear, feasible set is a polyhedron
opt_solvers()
#>     solver   package available                                 classes duals
#> 1    highs     highs      TRUE                              LP MILP QP  TRUE
#> 2 clarabel  clarabel      TRUE                        LP QP SOCP CONIC  TRUE
#> 3     scip      scip      TRUE LP MILP QP MIQP QCQP MIQCQP SOCP MISOCP FALSE
#> 4     osqp      osqp      TRUE                                   LP QP  TRUE
#> 5   nloptr    nloptr      TRUE                          NLP QCQP QP LP FALSE
#> 6      roi ROI, slam      TRUE                              LP MILP QP FALSE
#> 7    nlpbb    nloptr      TRUE                       MINLP MIQCQP MIQP FALSE
#>                                                                                    description
#> 1               HiGHS: dual simplex, interior point and branch and bound (LP, convex QP, MILP)
#> 2 Clarabel: interior point for LP, convex QP and SOC, PSD, exponential and power cone programs
#> 3         SCIP: branch and bound with spatial branching, global for nonconvex quadratic models
#> 4               OSQP: operator splitting for convex QP and LP (first-order, moderate accuracy)
#> 5                            NLopt SLSQP (default): local gradient-based nonlinear programming
#> 6               ROI: R Optimization Infrastructure with installed plugins (options$roi_solver)
#> 7   NLP branch and bound over NLopt relaxations (global for convex MINLP, heuristic otherwise)
```

[`opt_solve()`](https://adamrsophie.github.io/optcompR/reference/opt_solve.md)
with `solver = "auto"` picks the highest priority backend that supports
the problem class (here HiGHS for an LP).

``` r

r <- opt_solve(m)
r
#> <opt_result>
#>   solver:      highs 1.15.1.2
#>   problem:     LP (convex)
#>   status:      optimal (primal: feasible_point, dual: feasible_point)
#>   objective:   410
#>   time:        0.010s (build 0.001s)
opt_value(r, c(tables, chairs))
#> [1] 30 40
opt_objective_value(r)
#> [1] 410
```

## Duals, slacks and the solution table

``` r

opt_dual(r, carpentry)
#> carpentry 
#>       1.5
opt_dual(r, painting)
#> painting 
#>      0.5
opt_slack(r, carpentry)
#>         row value lower upper slack active
#> 1 carpentry     0  -Inf     0     0   TRUE
opt_solution(r)
#>     name value lower upper       type reduced_cost
#> 1 tables    30     0   Inf continuous            0
#> 2 chairs    40     0   Inf continuous            0
```

A dual is the derivative of the optimal objective with respect to the
bound of the row: one more carpentry hour raises profit by
`opt_dual(r, carpentry)`. This holds in both minimization and
maximization.

## Independent checks

[`opt_check()`](https://adamrsophie.github.io/optcompR/reference/opt_check.md)
recomputes row activities, bound and integrality violations and the KKT
residuals from the model itself, without trusting the solver.

``` r

opt_check(r)
#> <opt_diagnostics: feasible>
#>   max row violation (scaled):   0
#>   max bound violation (scaled): 0
#>   max integrality violation:    0
#>   objective recomputation error: 0
#>   KKT: stationarity 0, dual feasibility 0, complementarity 0 (pass)
```

## Where to go next

- [`vignette("linear-programming")`](https://adamrsophie.github.io/optcompR/articles/linear-programming.md),
  [`vignette("quadratic-programming")`](https://adamrsophie.github.io/optcompR/articles/quadratic-programming.md)
  and
  [`vignette("mixed-integer")`](https://adamrsophie.github.io/optcompR/articles/mixed-integer.md)
  for the main problem classes.
- [`vignette("parameterized-models")`](https://adamrsophie.github.io/optcompR/articles/parameterized-models.md)
  for repeated solves.
- [`vignette("nonlinear")`](https://adamrsophie.github.io/optcompR/articles/nonlinear.md)
  for nonlinear programming with exact derivatives.
- [`vignette("concepts")`](https://adamrsophie.github.io/optcompR/articles/concepts.md)
  for the mathematical background.
