# Constraint diagnostics

``` r

library(optcompR)
```

optcompR describes an optimization problem once, as variables,
expressions, constraints and an objective, and then classifies it,
chooses a capable solver Solvers work to tolerances, and their reported
status is a claim about the problem they received.
[`opt_check()`](https://adamrsophie.github.io/optcompR/reference/opt_check.md)
verifies a result against the original model.

## What is checked

- Row activities against their bounds, scaled by `max(1, |bound|)`.
- Variable bounds and integrality.
- The objective recomputed from the expression.
- KKT conditions for continuous problems: stationarity \nabla f - J^\top
  y - z = 0, sign conditions and complementarity.

Feasibility passes at `tol`; KKT conditions at `sqrt(tol)`.

``` r

m <- opt_model()
x <- opt_variable(m, "x", dim = 2, lower = 0)
opt_constraint(m, x[1] + 2 * x[2] >= 4, name = "cover")
opt_constraint(m, 3 * x[1] + x[2] >= 3, name = "mix")
opt_objective(m, x[1]^2 + x[2]^2)
r <- opt_solve(m)
d <- opt_check(r)
d
#> <opt_diagnostics: feasible>
#>   max row violation (scaled):   0
#>   max bound violation (scaled): 0
#>   max integrality violation:    0
#>   objective recomputation error: 0
#>   KKT: stationarity 5e-08, dual feasibility 0, complementarity 0 (pass)
d$rows
#>     row value lower upper violation
#> 1 cover     0     0   Inf         0
#> 2   mix     1     0   Inf         0
d$kkt
#> $stationarity
#> [1] 5e-08
#> 
#> $dual_feasibility
#> [1] 0
#> 
#> $complementarity
#> [1] 0
#> 
#> $stationarity_residual
#> [1] -8.0e-08 -1.6e-07
#> 
#> $pass
#> [1] TRUE
#> 
#> $estimated
#> [1] FALSE
```

## Results without multipliers

When a backend returns no duals (SCIP, NLopt),
[`opt_check()`](https://adamrsophie.github.io/optcompR/reference/opt_check.md)
estimates them by least squares on the active rows and bounds, and says
so:

``` r

r2 <- opt_solve(m, solver = "nloptr")
opt_check(r2)$kkt$estimated
#> [1] TRUE
```

## Slack and activity

``` r

opt_slack(r, "cover")
#>     row value lower upper slack active
#> 1 cover     0     0   Inf     0   TRUE
opt_slack(r, "mix")
#>   row value lower upper slack active
#> 1 mix     1     0   Inf     1  FALSE
```

## Derivatives at a point

The Jacobian of the constraint rows and the Hessian of the Lagrangian
are exact, so KKT residuals do not depend on finite differences:

``` r

opt_jacobian(m, r)
#> 2 x 2 sparse Matrix of class "dgCMatrix"
#>       x[1] x[2]
#> cover    1    2
#> mix      3    1
opt_hessian(m, r, multipliers = r$row_dual)
#> 2 x 2 sparse Matrix of class "dsCMatrix"
#>      x[1] x[2]
#> x[1]    2    .
#> x[2]    .    2
```
