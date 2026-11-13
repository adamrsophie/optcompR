# Robust and stochastic optimization

``` r

library(optcompR)
```

Uncertain data can be handled in two ways. Robust optimization asks a
constraint to hold for every value in an uncertainty set. Stochastic
optimization treats the data as random and constrains or optimizes
probabilities, expectations or risk measures. optcompR provides the
standard exact reformulations of both, built from ordinary variables,
rows and cones, so any capable solver handles the result. Each function
says whether its reformulation is exact or an approximation.

| Function | Model | Reformulation | Exact? |
|----|----|----|----|
| `opt_robust_constraint(..., "box")` | a \in \[a_0 - \delta, a_0 + \delta\] | linear | yes |
| `opt_robust_constraint(..., "budget")` | at most \Gamma deviations | linear (duality) | yes |
| `opt_robust_constraint(..., "ellipsoid")` | a = a_0 + Pu,\\ \\u\\ \le 1 | second-order cone | yes |
| `opt_chance_constraint(..., "gaussian")` | a \sim N(\mu, \Sigma) | second-order cone | yes, for p \ge 0.5 |
| `opt_chance_constraint(..., "scenario")` | sampled a | one row per sample | approximation |
| [`opt_cvar()`](https://adamrsophie.github.io/optcompR/reference/opt_cvar.md) | conditional value at risk | linear (Rockafellar and Uryasev) | yes |
| [`opt_expectation()`](https://adamrsophie.github.io/optcompR/reference/opt_cvar.md) | sample average | linear | yes for the sample |
| [`opt_worst_expectation()`](https://adamrsophie.github.io/optcompR/reference/opt_cvar.md) | finite ambiguity set | epigraph | yes |

## Robust linear constraints

A production plan uses two resources whose consumption per unit is
uncertain: nominally one unit each, but up to 0.5 and 0.2 more. With a
box uncertainty set the constraint must hold for the worst coefficients,
which for nonnegative production are the largest ones:

``` r

m <- opt_model("max")
x <- opt_variable(m, "x", dim = 2, lower = 0)
opt_robust_constraint(m, x, a0 = c(1, 1), b = 10, delta = c(0.5, 0.2))
opt_objective(m, 3 * x[1] + 2.8 * x[2])
r <- opt_solve(m)
opt_value(r, x)
#> [1] 0.000000 8.333333
opt_objective_value(r)
#> [1] 23.33333
```

The reformulation is a_0^\top x + \delta^\top s \le b with s \ge \|x\|,
so it is valid for variables of either sign.

### Budgeted uncertainty

Protecting against every coefficient at its worst at once is often too
conservative. The budget set of Bertsimas and Sim allows at most \Gamma
coefficients to deviate. \Gamma = 0 is the nominal problem and \Gamma =
n the box:

``` r

budget <- function(gamma) {
  m <- opt_model("max")
  x <- opt_variable(m, "x", dim = 3, lower = 0)
  opt_robust_constraint(m, x, a0 = c(1, 1, 1), b = 6, delta = c(1, 1, 1),
                        uncertainty = "budget", gamma = gamma)
  opt_objective(m, sum(x))
  opt_objective_value(opt_solve(m))
}
sapply(c(0, 0.5, 1, 2, 3), budget)
#> [1] 6.000000 5.142857 4.500000 3.600000 3.000000
```

### Ellipsoidal uncertainty

With a = a_0 + Pu and \\u\\\_2 \le 1, the worst case of a^\top x is
a_0^\top x + \\P^\top x\\\_2, a second-order cone:

``` r

m <- opt_model("max")
x <- opt_variable(m, "x", dim = 2, lower = 0)
P <- diag(c(0.5, 0.2))
opt_robust_constraint(m, x, a0 = c(1, 1), b = 10, uncertainty = "ellipsoid",
                      P = P)
opt_objective(m, 3 * x[1] + 2.8 * x[2])
r <- opt_solve(m)
v <- opt_value(r, x)
v
#> [1] 1.617836 6.799810
sum(v) + sqrt(sum((t(P) %*% v)^2))  # the worst case is exactly 10
#> [1] 10
opt_objective_value(r)
#> [1] 23.89298
```

The ellipsoid is smaller than the box with the same half-widths, so the
robust plan is less conservative.

## Chance constraints

A chance constraint asks P(a^\top x \le b) \ge p. For Gaussian
coefficients a \sim N(\mu, \Sigma) it is equivalent to \mu^\top x +
\Phi^{-1}(p) \\\Sigma^{1/2} x\\\_2 \le b, a convex second-order cone
when p \ge 0.5:

``` r

Sigma <- matrix(c(0.04, 0.01, 0.01, 0.02), 2)
m <- opt_model("max")
x <- opt_variable(m, "x", dim = 2, lower = 0)
opt_chance_constraint(m, x, b = 10, prob = 0.95, mean = c(1, 1),
                      Sigma = Sigma)
opt_objective(m, sum(x))
v <- opt_value(opt_solve(m), x)
v
#> [1] 2.053234 6.159687
pnorm((10 - sum(v)) / sqrt(drop(t(v) %*% Sigma %*% v)))
#> [1] 0.95
```

The achieved probability is exactly the requested 0.95. A Monte Carlo
check agrees:

``` r

set.seed(1)
a <- matrix(rnorm(2e5), ncol = 2) %*% chol(Sigma) + 1
mean(a %*% v <= 10)
#> [1] 0.94994
```

For other distributions the scenario approach imposes the constraint for
each of N samples. It is an approximation: with d decision variables the
solution satisfies the chance constraint at level 1 - \epsilon with
confidence 1 - \beta when N \ge \frac{2}{\epsilon}(\ln\frac1\beta + d)
(Calafiore and Campi). The bound is the user’s to check; optcompR does
not claim it.

``` r

set.seed(2)
S <- matrix(rexp(2 * 300), ncol = 2)  # skewed coefficients
m <- opt_model("max")
x <- opt_variable(m, "x", dim = 2, lower = 0)
opt_chance_constraint(m, x, b = 10, samples = S, method = "scenario")
opt_objective(m, sum(x))
v <- opt_value(opt_solve(m), x)
v
#> [1] 1.154559 1.542164
# out of sample probability of violation
mean(matrix(rexp(2e5), ncol = 2) %*% v > 10)
#> [1] 0.00558
```

## Risk measures over scenarios

### Sample average approximation

Given K loss scenarios L_k(x),
[`opt_expectation()`](https://adamrsophie.github.io/optcompR/reference/opt_cvar.md)
is their (weighted) average; minimizing it is sample average
approximation.

### Conditional value at risk

The conditional value at risk at level \alpha is the expected loss in
the worst 1 - \alpha fraction of outcomes. Rockafellar and Uryasev
showed

\mathrm{CVaR}\_\alpha(L) = \min_t\\ t + \frac{1}{1-\alpha} \sum_k p_k
(L_k - t)\_+ ,

which
[`opt_cvar()`](https://adamrsophie.github.io/optcompR/reference/opt_cvar.md)
writes with an auxiliary variable for t and one per scenario for the
excess. Minimizing it, or bounding it from above, is exact. A portfolio
over four return scenarios:

``` r

R <- matrix(c(0.10, 0.02, -0.05, 0.04,
              0.03, 0.01,  0.02, 0.01), 4,
            dimnames = list(NULL, c("stock", "bond")))
cvar_portfolio <- function(alpha) {
  m <- opt_model()
  w <- opt_variable(m, "w", dim = 2, lower = 0)
  opt_constraint(m, sum(w) == 1)
  loss <- -(R %*% w)
  opt_objective(m, opt_cvar(m, loss, alpha = alpha))
  r <- opt_solve(m)
  c(opt_value(r, w), cvar = opt_objective_value(r))
}
rbind(`alpha = 0.25` = cvar_portfolio(0.25),
      `alpha = 0.75` = cvar_portfolio(0.75))
#>                                 cvar
#> alpha = 0.25 0.000 1.000 -0.01333333
#> alpha = 0.75 0.125 0.875 -0.01125000
```

At \alpha = 0.75 only the worst scenario counts and the portfolio moves
toward bonds. CVaR can also be a constraint: maximize expected return
subject to a CVaR limit.

``` r

m <- opt_model("max")
w <- opt_variable(m, "w", dim = 2, lower = 0)
opt_constraint(m, sum(w) == 1)
opt_constraint(m, opt_cvar(m, -(R %*% w), alpha = 0.75) <= 0.01)
opt_objective(m, opt_expectation(R %*% w))
opt_value(opt_solve(m), w)
#> [1] 0.4285714 0.5714286
```

A CVaR expression must appear where a larger value is worse (minimized,
or on the left of `<=`). Maximizing it is not modeled by this
formulation, because the auxiliary variables can then inflate it.

### Distributionally robust expectations

When the scenario probabilities are themselves uncertain,
[`opt_worst_expectation()`](https://adamrsophie.github.io/optcompR/reference/opt_cvar.md)
takes a matrix of candidate distributions and returns the largest
expected loss among them, which is also the worst case over their convex
hull:

``` r

probs <- rbind(uniform = rep(0.25, 4),
               crash   = c(0.1, 0.1, 0.7, 0.1),
               boom    = c(0.7, 0.1, 0.1, 0.1))
m <- opt_model()
w <- opt_variable(m, "w", dim = 2, lower = 0)
opt_constraint(m, sum(w) == 1)
opt_objective(m, opt_worst_expectation(m, -(R %*% w), probs))
r <- opt_solve(m)
opt_value(r, w)
#> [1] 0.03125 0.96875
probs %*% (-(R %*% opt_value(r, w)))
#>               [,1]
#> uniform -0.0178125
#> crash   -0.0178125
#> boom    -0.0264375
```

The optimal weights equalize the expected loss under the two binding
distributions.
