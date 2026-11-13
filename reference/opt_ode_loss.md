# Least squares loss of an ODE model from modelcompR

Builds the scalar expression \\\sum\_{t, s} w\_{ts} (u_s(t; \theta) -
y\_{ts})^2\\ for an `ODEProblem` of the modelcompR package, with
\\\theta\\ given by decision variables. Each evaluation solves the ODE
and its forward sensitivities (`ODEForwardSensitivityProblem()`), so the
gradient is the exact derivative of the numerical solution up to the
integration tolerances; no finite differences are used. Value and
gradient at the same point share one solve.

## Usage

``` r
opt_ode_loss(
  prob,
  theta,
  times,
  data,
  states = NULL,
  weights = NULL,
  reltol = 1e-10,
  abstol = 1e-12,
  alg = NULL
)
```

## Arguments

- prob:

  A modelcompR `ODEProblem`; its parameters are replaced by `theta`.

- theta:

  Decision variables (or an expression) for the ODE parameters, in the
  order the right-hand side uses them.

- times:

  Observation times.

- data:

  Matrix of observations, one row per time and one column per observed
  state (`NA` for missing values).

- states:

  Indices of the observed states.

- weights:

  Matrix of weights like `data` (default 1).

- reltol, abstol:

  Integration tolerances.

- alg:

  Integration algorithm (default
  [`modelcompR::Tsit5()`](https://rdrr.io/pkg/modelcompR/man/ode_algorithms.html)).

## Value

A scalar nonlinear `opt_expr`.

## Details

The loss is generally not convex, so solves are local: start from a
sensible `start` on the variables.

## Examples

``` r
f <- function(u, p, t) -p[1] * u
prob <- modelcompR::ODEProblem(f, 1, c(0, 2), 1)
times <- c(0.5, 1, 1.5, 2)
y <- matrix(exp(-0.8 * times), ncol = 1)
m <- opt_model()
k <- opt_variable(m, "k", lower = 0, upper = 5, start = 0.3)
opt_objective(m, opt_ode_loss(prob, k, times, y))
opt_value(opt_solve(m), k)  # 0.8
#> [1] 0.8
```
