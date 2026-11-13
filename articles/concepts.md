# Optimization concepts

``` r

library(optcompR)
```

optcompR describes an optimization problem once, as variables,
expressions, constraints and an objective, and then classifies it,
chooses a capable solver This tutorial explains the distinctions that
decide which solver applies and how to read its result, with a small
model for each.

## Convex and nonconvex

A problem is convex when it minimizes a convex function (or maximizes a
concave one) over a convex set. Every local minimum of a convex problem
is global.
[`opt_classify()`](https://adamrsophie.github.io/optcompR/reference/opt_classify.md)
decides convexity for polynomial models of degree two and records how:

``` r

m <- opt_model()
x <- opt_variable(m, "x", dim = 2)
opt_objective(m, x[1]^2 + x[1] * x[2] + x[2]^2)
cl <- opt_classify(m)
cl$convexity
#> [1] "convex"
cl$evidence
#> [1] "syntactic: objective is quadratic; constraint rows: none"                            
#> [2] "numerical eigenvalues: objective is convex (smallest eigenvalue 1 (tolerance 3e-09))"
opt_objective(m, x[1]^2 + 3 * x[1] * x[2] + x[2]^2)
opt_classify(m)$convexity
#> [1] "nonconvex"
```

For general nonlinear expressions convexity is reported as `unknown`; it
is not guessed. Cone constraints are convex sets by construction, so a
model whose only nonlinear structure is cones (with a linear or convex
quadratic objective) is classified as convex `"SOCP"` or `"CONIC"`; see
[`vignette("conic")`](https://adamrsophie.github.io/optcompR/articles/conic.md).

## Continuous and mixed-integer

Integer variables make the feasible set a union of points or slices. The
LP relaxation gives a bound, and branch and bound closes the gap:

``` r

m <- opt_model("max")
x <- opt_variable(m, "x", dim = 2, type = "integer", lower = 0)
opt_constraint(m, 2 * x[1] + 2 * x[2] <= 7)
opt_objective(m, x[1] + x[2])
r_int <- opt_solve(m)
opt_set_type(m, x, "continuous")
r_lp <- opt_solve(m)
c(integer = opt_objective_value(r_int), relaxation = opt_objective_value(r_lp))
#>    integer relaxation 
#>        3.0        3.5
```

## Local and global optima

``` r

m <- opt_model()
x <- opt_variable(m, "x", lower = -3, upper = 3, start = 2)
opt_objective(m, x^4 - 4 * x^2 + x)
r <- opt_solve(m, solver = "nloptr")
c(x = opt_value(r, x), f = opt_objective_value(r))
#>         x         f 
#>  1.346997 -2.618556
opt_status(r)
#> [1] "locally_optimal"
opt_set_start(m, x, -2)
r2 <- opt_solve(m, solver = "nloptr")
c(x = opt_value(r2, x), f = opt_objective_value(r2))
#>         x         f 
#> -1.472998 -5.444192
```

Different starts end in different local minima, and the status says
`locally_optimal`, never `optimal`.

## Feasible and optimal

A feasible point satisfies every constraint; an optimal point is
feasible and has the best objective. Solvers stopped by limits may
return a feasible point that is not proven optimal;
[`opt_has_solution()`](https://adamrsophie.github.io/optcompR/reference/opt_status.md),
[`opt_status()`](https://adamrsophie.github.io/optcompR/reference/opt_status.md)
and
[`opt_gap()`](https://adamrsophie.github.io/optcompR/reference/opt_objective_value.md)
together tell which case applies.

## Primal and dual

The primal variables are the decisions x. The dual variables y price the
constraints: y_i = \partial f^\* / \partial b_i. In optcompR this sign
convention holds for both minimization and maximization, and
[`opt_check()`](https://adamrsophie.github.io/optcompR/reference/opt_check.md)
verifies the KKT stationarity condition \nabla f(x) - J(x)^\top y - z =
0.

## Termination and solution status

[`opt_status()`](https://adamrsophie.github.io/optcompR/reference/opt_status.md)
reports why the solver stopped (`optimal`, `infeasible`, `time_limit`,
…).
[`opt_primal_status()`](https://adamrsophie.github.io/optcompR/reference/opt_status.md)
and
[`opt_dual_status()`](https://adamrsophie.github.io/optcompR/reference/opt_status.md)
report whether the returned vectors are feasible. A `time_limit`
termination with a `feasible_point` primal status is a usable incumbent.

## Derivatives and sensitivity of optima

The gradient of the objective, \nabla f(x), describes the model
function. The sensitivity of the optimum, \partial f^\* / \partial b,
describes the optimization problem as a whole; by the envelope theorem
it equals the dual. See
[`vignette("sensitivity")`](https://adamrsophie.github.io/optcompR/articles/sensitivity.md).
The derivative of the optimal point itself, \partial x^\* / \partial p,
follows from differentiating the KKT conditions and is computed by
[`opt_solution_jacobian()`](https://adamrsophie.github.io/optcompR/reference/opt_solution_jacobian.md);
see
[`vignette("solution-derivatives")`](https://adamrsophie.github.io/optcompR/articles/solution-derivatives.md).
