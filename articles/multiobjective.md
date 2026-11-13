# Multiobjective optimization

``` r

library(optcompR)
```

optcompR describes an optimization problem once, as variables,
expressions, constraints and an objective, and then classifies it,
chooses a capable solver With several objectives there is usually no
single best point, only a set of Pareto optimal points where no
objective can improve without another getting worse. optcompR offers
three scalarizations. Each restores the model’s own objective and
removes temporary constraints afterwards.

``` r

m <- opt_model()
x <- opt_variable(m, "x", lower = 0, upper = 4)
y <- opt_variable(m, "y", lower = 0, upper = 4)
opt_constraint(m, x + y >= 4)
cost <- 2 * x + y
risk <- x + 3 * y
```

## Weighted sums

``` r

r <- opt_weighted(m, list(cost = cost, risk = risk), weights = c(1, 1))
r$objective_values
#> cost risk 
#>    8    4
```

Weighted sums with positive weights give Pareto optimal points, but only
supported ones: points on nonconvex parts of the front (and many points
of integer problems) are never reached.

## Lexicographic order

``` r

st <- opt_lexicographic(m, list(cost = cost, risk = risk))
st[[2]]$objective_values
#> cost risk 
#>    4   12
opt_value(st[[2]], c(x, y))
#> [1] 4e-07 4e+00
```

Cost is optimized first; risk is then minimized among points whose cost
is within a relative tolerance `tol` of the best.

## Epsilon constraints and the Pareto front

``` r

front <- opt_epsilon_constraint(m, cost, risk, epsilon = seq(4, 12, by = 2))
front
#>   epsilon  status primary secondary
#> 1       4 optimal       8         4
#> 2       6 optimal       7         6
#> 3       8 optimal       6         8
#> 4      10 optimal       5        10
#> 5      12 optimal       4        12
plot(front$secondary, front$primary, type = "b",
     xlab = "risk", ylab = "cost", main = "Pareto front")
```

![](multiobjective_files/figure-html/unnamed-chunk-5-1.png)

Epsilon constraints reach every Pareto point, including unsupported
ones, but can return weakly Pareto optimal points when the bound on the
secondary objective is not tight.
