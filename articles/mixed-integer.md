# Mixed-integer programming

``` r

library(optcompR)
```

optcompR describes an optimization problem once, as variables,
expressions, constraints and an objective, and then classifies it,
chooses a capable solver Integer and binary variables turn an LP into a
MILP (or a QP into a MIQP). Branch and bound solvers report an incumbent
objective, a proven bound and the relative gap between them.

## Knapsack

``` r

value <- c(10, 13, 7, 8, 12, 4)
weight <- c(5, 7, 3, 4, 6, 2)
m <- opt_model("max")
take <- opt_variable(m, "take", dim = 6, type = "binary")
cap <- opt_constraint(m, sum(weight * take) <= 15, name = "capacity")
opt_objective(m, sum(value * take))
r <- opt_solve(m)
opt_status(r)
#> [1] "optimal"
which(opt_value(r, take) > 0.5)
#> [1] 3 4 5 6
c(objective = opt_objective_value(r), bound = opt_objective_bound(r),
  gap = opt_gap(r))
#> objective     bound       gap 
#>        31        31         0
```

## Facility location with fixed costs

Binary open decisions link to continuous flows through big M
constraints.

``` r

fixed <- c(30, 25, 40)
unit <- matrix(c(2, 4, 5,
                 3, 1, 3,
                 4, 3, 1,
                 5, 4, 2), 4, byrow = TRUE)  # 4 customers by 3 sites
demand <- c(10, 15, 12, 8)
m2 <- opt_model()
open <- opt_variable(m2, "open", dim = 3, type = "binary")
flow <- opt_variable(m2, "flow", dim = c(4, 3), lower = 0)
for (i in 1:4) opt_constraint(m2, sum(flow[i, ]) == demand[i])
for (j in 1:3) opt_constraint(m2, sum(flow[, j]) <= sum(demand) * open[j])
opt_objective(m2, sum(fixed * open) + sum(unit * flow))
r2 <- opt_solve(m2)
opt_value(r2, open)
#> [1] 0 1 1
opt_objective_value(r2)
#> [1] 148
```

## Limits and gaps

``` r

r3 <- opt_solve(m2, mip_gap = 0.05, time_limit = 10)
opt_status(r3)
#> [1] "optimal"
opt_gap(r3)
#> [1] 0
```

A `time_limit` status with a feasible point is a usable but unproven
solution:
[`opt_has_solution()`](https://adamrsophie.github.io/optcompR/reference/opt_status.md)
tells whether a point exists, and the gap tells how far it may be from
optimal.

## Duals of a MILP

Integer programs have no duals in the LP sense.
[`opt_fixed_relaxation()`](https://adamrsophie.github.io/optcompR/reference/opt_fixed_relaxation.md)
fixes the integer variables at their optimal values and solves the
remaining LP, whose duals are the usual marginal values for that fixed
configuration.

``` r

fr <- opt_fixed_relaxation(r)
opt_dual(fr, cap)
#> capacity 
#>        0
```

## Mixed-integer quadratic problems

``` r

m4 <- opt_model()
x <- opt_variable(m4, "x", dim = 3, type = "integer", lower = -5, upper = 5)
opt_objective(m4, sum((x - c(1.4, -2.6, 0.5))^2))
r4 <- opt_solve(m4)
r4$solver
#> [1] "scip"
opt_value(r4, x)
#> [1]  1 -3  0
```
