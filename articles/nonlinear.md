# Nonlinear optimization

``` r

library(optcompR)
```

optcompR describes an optimization problem once, as variables,
expressions, constraints and an objective, and then classifies it,
chooses a capable solver Nonlinear expressions combine variables with
`exp`, `log`, `sqrt`, `sin`, `cos`, `tan`, `tanh`, `sinh`, `cosh`,
`abs`, `log1p`, `expm1`, `atan`, `asin`, `acos`, powers, products and
quotients. Every expression has an exact gradient, Jacobian and Hessian
built by the chain rule; there are no finite differences unless you ask
for them.

## Rosenbrock

``` r

m <- opt_model()
x <- opt_variable(m, "x", dim = 2, start = c(-1.2, 1))
opt_objective(m, (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2)
opt_classify(m)
#> <opt_classification: NLP, unknown>
#>   - syntactic: objective is nonlinear; constraint rows: none
#>   - unknown: convexity of general nonlinear expressions is not analyzed; solutions from local solvers are local optima
r <- opt_solve(m)
opt_status(r)
#> [1] "locally_optimal"
opt_value(r, x)
#> [1] 1 1
```

The status is `locally_optimal`: SLSQP is a local method, and the
problem is not known to be convex.

## Hock and Schittkowski problem 71

\min\\ x_1 x_4 (x_1 + x_2 + x_3) + x_3 subject to x_1 x_2 x_3 x_4 \ge
25, \sum x_i^2 = 40 and 1 \le x_i \le 5.

``` r

m <- opt_model()
x <- opt_variable(m, "x", dim = 4, lower = 1, upper = 5, start = c(1, 5, 5, 1))
opt_constraint(m, x[1] * x[2] * x[3] * x[4] >= 25, name = "product")
opt_constraint(m, sum(x^2) == 40, name = "sphere")
opt_objective(m, x[1] * x[4] * (x[1] + x[2] + x[3]) + x[3])
r <- opt_solve(m)
opt_objective_value(r)  # reference 17.0140173
#> [1] 17.01402
opt_check(r)$kkt
#> $stationarity
#> [1] 4.782548e-09
#> 
#> $dual_feasibility
#> [1] 0
#> 
#> $complementarity
#> [1] 1.405089e-11
#> 
#> $stationarity_residual
#> [1]  2.220446e-16  6.597087e-09 -8.619240e-09  1.192852e-09
#> 
#> $pass
#> [1] TRUE
#> 
#> $estimated
#> [1] TRUE
```

## IPOPT

[IPOPT](https://coin-or.github.io/Ipopt/) is an interior point solver
for large nonlinear programs. optcompR binds to its C interface when
IPOPT is found with `pkg-config` at installation (`brew install ipopt`
on macOS, `apt install coinor-libipopt-dev` on Debian and Ubuntu), and
then
[`opt_ipopt_available()`](https://adamrsophie.github.io/optcompR/reference/opt_ipopt_available.md)
is `TRUE`. IPOPT receives the exact gradient, Jacobian and Hessian of
the Lagrangian; for external functions without a Hessian it uses a
limited-memory quasi-Newton approximation. Unlike SLSQP it returns
multipliers, so duals and the KKT check use solver values rather than
estimates.

``` r

r <- opt_solve(m, solver = "ipopt")
opt_status(r)
#> [1] "locally_optimal"
opt_objective_value(r)
#> [1] 17.01402
opt_dual(r, "product")
#>   product 
#> 0.5522937
opt_check(r)$kkt$estimated
#> [1] FALSE
```

IPOPT is a local solver: on nonconvex models its result is
`"locally_optimal"`, exactly as for NLopt. The binding passes dense
Jacobian and Hessian structures, which suits models with up to a few
hundred variables.

## Exact derivatives

``` r

m <- opt_model()
x <- opt_variable(m, "x")
y <- opt_variable(m, "y")
opt_objective(m, x^2 + 3 * x * y + exp(y))
opt_constraint(m, x * y <= 4)
opt_gradient(m, c(1, 2))
#>        x        y 
#>  8.00000 10.38906
opt_jacobian(m, c(1, 2))
#> 1 x 2 sparse Matrix of class "dgCMatrix"
#>    x y
#> c1 2 1
opt_hessian(m, c(1, 2))
#> 2 x 2 sparse Matrix of class "dsCMatrix"
#>   x        y
#> x 2 3.000000
#> y 3 7.389056
```

## External R functions

Functions that are not written with optcompR operators, such as
simulation losses, enter through
[`opt_function()`](https://adamrsophie.github.io/optcompR/reference/opt_function.md).
The gradient source is explicit:

``` r

rosen <- function(v) (1 - v[1])^2 + 100 * (v[2] - v[1]^2)^2
rosen_grad <- function(v) c(-2 * (1 - v[1]) - 400 * v[1] * (v[2] - v[1]^2),
                            200 * (v[2] - v[1]^2))
m <- opt_model()
x <- opt_variable(m, "x", dim = 2, start = c(-1.2, 1))
opt_objective(m, opt_function(rosen, x, gradient = rosen_grad, name = "rosen"))
opt_value(opt_solve(m), x)
#> [1] 1 1
```

With `gradient = "finite_difference"` the gradient is approximated by
central differences; results then carry a warning and list the function
in `r$metadata$finite_difference`. With `gradient = "jaxR"` and the jaxR
package, the gradient comes from reverse-mode automatic differentiation
of `fn`.

A model from modelcompR (for example the sum of squared errors of an ODE
solution against data, as a function of its parameters) can be optimized
the same way.
[`opt_ode_loss()`](https://adamrsophie.github.io/optcompR/reference/opt_ode_loss.md)
builds that loss with gradients from modelcompR’s forward sensitivities;
see
[`vignette("ode-fitting")`](https://adamrsophie.github.io/optcompR/articles/ode-fitting.md).

## Limits

- NLopt SLSQP finds local solutions.
- Mixed-integer nonlinear problems (MINLP) are solved by NLP branch and
  bound (`solver = "nlpbb"`); results are global only when the
  continuous relaxation is known to be convex, and labeled local
  otherwise. See
  [`vignette("discrete-structure")`](https://adamrsophie.github.io/optcompR/articles/discrete-structure.md).
- Convex problems with norms, logarithms, exponentials or matrix
  eigenvalues are often better written as conic programs, which Clarabel
  solves to global optimality; see
  [`vignette("conic")`](https://adamrsophie.github.io/optcompR/articles/conic.md).
