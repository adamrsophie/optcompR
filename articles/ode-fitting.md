# Fitting ODE models with modelcompR

``` r

library(optcompR)
```

[modelcompR](https://github.com/adamrsophie/modelcompR) solves
differential equations in R and computes their sensitivities with
respect to parameters.
[`opt_ode_loss()`](https://adamrsophie.github.io/optcompR/reference/opt_ode_loss.md)
turns an ODE model into a least squares objective whose gradient comes
from forward sensitivities, so a gradient-based solver can fit ODE
parameters to data with exact derivatives of the numerical solution.

## A two-compartment model

Drug in a depot compartment u_1 is absorbed at rate k_a into a central
compartment u_2, which is eliminated at rate k_e:

u_1' = -k_a u_1, \qquad u_2' = k_a u_1 - k_e u_2 .

``` r

f <- function(u, p, t) c(-p[1] * u[1], p[1] * u[1] - p[2] * u[2])
prob <- modelcompR::ODEProblem(f, c(1, 0), c(0, 6), c(0.5, 0.5))
```

Simulated observations of both compartments, with noise:

``` r

truth <- c(ka = 0.9, ke = 0.25)
times <- c(0.25, 0.5, 1, 1.5, 2, 3, 4, 6)
exact <- function(th, t) {
  u1 <- exp(-th[1] * t)
  u2 <- th[1] / (th[2] - th[1]) * (exp(-th[1] * t) - exp(-th[2] * t))
  cbind(u1, u2)
}
set.seed(4)
y <- exact(truth, times) + matrix(rnorm(16, sd = 0.01), 8)
```

## The fit

The parameters are decision variables with bounds and a starting point.
The loss is \sum\_{t, s} w\_{ts} (u_s(t; \theta) - y\_{ts})^2.

``` r

m <- opt_model()
theta <- opt_variable(m, "theta", dim = 2, lower = 0.01, upper = 5,
                      start = c(0.3, 0.6))
opt_objective(m, opt_ode_loss(prob, theta, times, y))
opt_classify(m)$class
#> [1] "NLP"
r <- opt_solve(m)
opt_status(r)
#> [1] "locally_optimal"
opt_value(r, theta)
#> [1] 0.9011146 0.2485141
```

The status is `"locally_optimal"`: least squares fits of ODE models are
not convex in general, and the solver reports a local optimum. Starting
from several points is a cheap check.

## The gradient

Each evaluation of the loss solves the ODE together with its forward
sensitivity equations. The gradient

\nabla\_\theta \ell = 2 \sum\_{t, s} w\_{ts} (u_s(t) - y\_{ts})
\frac{\partial u_s(t)}{\partial \theta}

is assembled from those sensitivities. It agrees with finite differences
to the integration tolerances:

``` r

at <- c(0.6, 0.4)
g <- opt_gradient(m, at)
loss <- function(v) optcompR:::eval_expr(m$objective, v)
num <- vapply(1:2, function(j) {
  h <- 1e-6
  e <- replace(c(0, 0), j, h)
  (loss(at + e) - loss(at - e)) / (2 * h)
}, 0)
rbind(sensitivities = unname(g), finite_differences = num)
#>                         [,1]      [,2]
#> sensitivities      -1.219944 0.9196711
#> finite_differences -1.219944 0.9196712
```

## Partial observations and weights

Observe only the central compartment by passing `states`, give each
observation a weight, and mark missing values with `NA`:

``` r

y2 <- y[, 2, drop = FALSE]
y2[3] <- NA
m2 <- opt_model()
theta2 <- opt_variable(m2, "theta", dim = 2, lower = 0.01, upper = 5,
                       start = c(0.3, 0.6))
opt_objective(m2, opt_ode_loss(prob, theta2, times, y2, states = 2,
                               weights = matrix(1 / 0.01^2, 8, 1)))
opt_value(opt_solve(m2), theta2)
#> [1] 0.9387218 0.2498122
```

## Constraints and priors

Because the loss is an ordinary optcompR expression, the fit can carry
constraints and penalty terms. A ridge penalty toward prior values and a
constraint that absorption is faster than elimination:

``` r

m3 <- opt_model()
th <- opt_variable(m3, "theta", dim = 2, lower = 0.01, upper = 5,
                   start = c(0.6, 0.3))
opt_constraint(m3, th[1] - th[2] >= 0.1)
opt_objective(m3, opt_ode_loss(prob, th, times, y) +
                0.001 * sum((th - c(1, 0.2))^2))
opt_value(opt_solve(m3), th)
#> [1] 0.9012141 0.2485104
```
