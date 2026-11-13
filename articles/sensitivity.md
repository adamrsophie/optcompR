# Sensitivity and derivatives

``` r

library(optcompR)
```

optcompR describes an optimization problem once, as variables,
expressions, constraints and an objective, and then classifies it,
chooses a capable solver Two different questions are often both called
sensitivity:

1.  How do the model functions change with the decision variables? This
    is answered by derivatives of f and g
    ([`opt_gradient()`](https://adamrsophie.github.io/optcompR/reference/opt_gradient.md),
    [`opt_jacobian()`](https://adamrsophie.github.io/optcompR/reference/opt_gradient.md),
    [`opt_hessian()`](https://adamrsophie.github.io/optcompR/reference/opt_gradient.md)).
2.  How does the optimal value or solution change with the data? This is
    answered by duals, ranging and re-solving.

## Duals as derivatives of the optimal value

``` r

m <- opt_model("max")
x <- opt_variable(m, "x", lower = 0)
y <- opt_variable(m, "y", lower = 0)
cap <- opt_constraint(m, x + y <= 4, name = "cap")
lim <- opt_constraint(m, x <= 2, name = "lim")
opt_objective(m, 3 * x + 2 * y)
r <- opt_solve(m)
opt_dual(r, cap)
#> cap 
#>   2
```

The dual equals the finite difference of the optimal value:

``` r

b <- opt_parameter(m, "b", value = 4)
opt_delete_constraint(m, "cap")
opt_constraint(m, x + y <= b, name = "cap")
f <- opt_scenarios(m, b, c(4, 4.01))
diff(f$objective) / 0.01
#> [1] 2
```

## Ranging

For LPs solved by HiGHS, ranging gives the interval over which each cost
coefficient or row bound can move with the same optimal basis.

``` r

s <- opt_sensitivity(opt_solve(m))
s$costs
#>   variable coefficient lower upper objective_at_lower objective_at_upper
#> 1        x           3     2   Inf                  8                Inf
#> 2        y           2     0     3                  6                 12
s$rows
#>   row row_lower row_upper dual bound_lower bound_upper objective_at_lower
#> 1 lim      -Inf         2    1           0           4                  8
#> 2 cap      -Inf         4    2           2         Inf                  6
#>   objective_at_upper
#> 1                 12
#> 2                Inf
```

## Scenarios beyond the ranging interval

Outside those intervals the optimal value is still piecewise linear in a
right hand side, but with a different slope:

``` r

opt_scenarios(m, b, seq(1, 7, by = 1))
#>   scenario value  status objective
#> 1        1     1 optimal         3
#> 2        2     2 optimal         6
#> 3        3     3 optimal         8
#> 4        4     4 optimal        10
#> 5        5     5 optimal        12
#> 6        6     6 optimal        14
#> 7        7     7 optimal        16
```

## Mixed-integer models

``` r

mi <- opt_model("max")
z <- opt_variable(mi, "z", type = "integer", lower = 0, upper = 3)
w <- opt_variable(mi, "w", lower = 0)
cp <- opt_constraint(mi, z + w <= 4.5, name = "cap")
opt_objective(mi, 3 * z + 2 * w)
ri <- opt_solve(mi)
opt_dual(opt_fixed_relaxation(ri), cp)
#> cap 
#>   2
```

The duals of the fixed relaxation are valid for small changes that keep
the integer part of the solution optimal.

## Derivatives of model functions

``` r

q <- opt_model()
u <- opt_variable(q, "u", dim = 2)
opt_objective(q, u[1]^2 * u[2] + sin(u[2]))
opt_gradient(q, c(1, 0.5))
#>     u[1]     u[2] 
#> 1.000000 1.877583
opt_hessian(q, c(1, 0.5))
#> 2 x 2 sparse Matrix of class "dsCMatrix"
#>      u[1]       u[2]
#> u[1]    1  2.0000000
#> u[2]    2 -0.4794255
```

Derivatives of the optimal solution with respect to parameters
(differentiable optimization) are computed by
[`opt_solution_jacobian()`](https://adamrsophie.github.io/optcompR/reference/opt_solution_jacobian.md);
see
[`vignette("solution-derivatives")`](https://adamrsophie.github.io/optcompR/articles/solution-derivatives.md).
