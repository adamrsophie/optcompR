# Quadratic programming

``` r

library(optcompR)
```

A quadratic program has a quadratic objective and linear constraints:

\min_x\\ \tfrac12 x^\top Q x + c^\top x \quad \text{subject to}\quad l
\le Ax \le u.

It is convex when Q is positive semidefinite. optcompR stores Q in this
one-half convention and checks convexity before choosing a solver.

## Portfolio selection

``` r

mu <- c(0.08, 0.12, 0.10)
Sigma <- matrix(c(0.10, 0.02, 0.04,
                  0.02, 0.15, 0.03,
                  0.04, 0.03, 0.12), 3)
m <- opt_model()
w <- opt_variable(m, "w", dim = 3, lower = 0)
opt_constraint(m, sum(w) == 1, name = "budget")
ret <- opt_constraint(m, sum(mu * w) >= 0.10, name = "return")
opt_objective(m, t(w) %*% Sigma %*% w)
opt_classify(m)
#> <opt_classification: QP, convex>
#>   - syntactic: objective is quadratic; constraint rows: linear = 2
#>   - numerical eigenvalues: objective is convex (smallest eigenvalue 0.137316 (rounding level 2.48609e-14, tolerance 3.73211e-10))
r <- opt_solve(m)
round(opt_value(r, w), 4)
#> [1] 0.3469 0.3469 0.3061
opt_objective_value(r)
#> [1] 0.06102041
opt_dual(r, ret)
#>    return 
#> 0.7142857
```

The dual of the return constraint is the marginal variance of asking for
more return.

## Convexity evidence

``` r

opt_classify(m)$evidence
#> [1] "syntactic: objective is quadratic; constraint rows: linear = 2"                                                               
#> [2] "numerical eigenvalues: objective is convex (smallest eigenvalue 0.137316 (rounding level 2.48609e-14, tolerance 3.73211e-10))"
```

## A nonconvex QP

If Q is indefinite, local solvers can stop at a local minimum. optcompR
then prefers a global solver (SCIP) and refuses solvers that only accept
convex QPs.

``` r

m2 <- opt_model()
x <- opt_variable(m2, "x", dim = 2, lower = -1, upper = 1)
opt_objective(m2, x[1]^2 - x[2]^2 + 0.5 * x[1] * x[2])
opt_classify(m2)$convexity
#> [1] "nonconvex"
r2 <- opt_solve(m2)
r2$solver
#> [1] "scip"
opt_value(r2, x)
#> [1] -0.2490238  1.0000000
```

## Solver agreement

``` r

r_highs <- opt_solve(m, solver = "highs")
r_osqp <- opt_solve(m, solver = "osqp")
rbind(highs = opt_value(r_highs, w), osqp = opt_value(r_osqp, w))
#>            [,1]      [,2]      [,3]
#> highs 0.3469388 0.3469388 0.3061225
#> osqp  0.3469388 0.3469388 0.3061224
```

## Exact derivatives

``` r

opt_gradient(m, r)
#>      w[1]      w[2]      w[3] 
#> 0.1077551 0.1363265 0.1220408
opt_hessian(m, r)
#> 3 x 3 sparse Matrix of class "dsCMatrix"
#>      w[1] w[2] w[3]
#> w[1] 0.20 0.04 0.08
#> w[2] 0.04 0.30 0.06
#> w[3] 0.08 0.06 0.24
```

The Hessian equals Q (here 2\Sigma, because the objective is w^\top
\Sigma w).
