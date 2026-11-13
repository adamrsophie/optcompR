# optcompR <img src="man/figures/logo.png" align="right" height="139" alt="optcompR logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/adamrsophie/optcompR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/adamrsophie/optcompR/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/adamrsophie/optcompR/actions/workflows/pkgdown.yaml/badge.svg)](https://adamrsophie.github.io/optcompR/)
<!-- badges: end -->

Unified mathematical optimization and computational modeling in R.

optcompR is a solver independent modeling layer. Variables, parameters,
constraints and objectives are written with ordinary R operators; models are
classified (LP, MILP, QP, MIQP, QCQP, MIQCQP, SOCP, CONIC, NLP, MINLP) with
recorded convexity
evidence, lowered to sparse standard forms and solved by installed backends.
Results carry standardized statuses, duals with one documented sign
convention, and diagnostics recomputed independently from the model.

```r
library(optcompR)

m <- opt_model("max")
x <- opt_variable(m, "x", lower = 0)
y <- opt_variable(m, "y", lower = 0)
cap <- opt_constraint(m, x + y <= 4, name = "capacity")
opt_constraint(m, x <= 2)
opt_constraint(m, y <= 3)
opt_objective(m, 3 * x + 2 * y)

r <- opt_solve(m)
opt_status(r)          # "optimal"
opt_value(r, c(x, y))  # 2 2
opt_dual(r, cap)       # shadow price 2
opt_check(r)           # independent feasibility and KKT checks
```

## Features

* Vectorized expressions with sparse `Matrix` support, built in milliseconds
  for 1e5 terms.
* Parameters for reusable models, scenario sweeps and warm starts.
* Exact gradients, Jacobians and Lagrangian Hessians of every model function.
* External R functions with supplied, jaxR or explicitly labeled finite
  difference gradients.
* Sensitivity: duals as derivatives of the optimal value, HiGHS ranging,
  fixed relaxations of MILPs.
* Multiobjective scalarizations: weighted sums, lexicographic order and
  epsilon constraints.
* Exact feasibility and optimality certificates for LP and MILP through Z3
  (z3r).
* Second-order, semidefinite, exponential and power cone constraints.
* Indicator constraints, SOS1 and SOS2 sets, complementarity and logical
  constraints, native in SCIP and exactly reformulated elsewhere.
* Mixed-integer nonlinear programming by NLP branch and bound.
* Robust constraints (box, budget, ellipsoid), Gaussian and scenario chance
  constraints, CVaR, sample average and distributionally robust
  expectations.
* Derivatives of optimal solutions with respect to parameters by implicit
  differentiation of the KKT conditions.
* Least squares fitting of modelcompR ODE models with exact sensitivities.
* Conversion to and from symcompR expressions.

## Solver support

| Backend | Package | Classes | Global for | Duals | Ranging | Warm start |
|---|---|---|---|---|---|---|
| highs | highs | LP, MILP, convex QP | all of these | LP, QP | LP | yes |
| clarabel | clarabel | LP, convex QP, SOCP, CONIC (SOC, PSD, exp, power) | all of these | yes | no | no |
| scip | scip | LP, MILP, QP, MIQP, QCQP, MIQCQP, SOCP, MISOCP | all, including nonconvex | no | no | no |
| osqp | osqp | convex QP, LP | convex problems | yes | no | no |
| nloptr | nloptr | NLP, QCQP, QP, LP (SLSQP) | convex problems | estimated | no | start point |
| ipopt | IPOPT library | NLP, QCQP, QP, LP (exact Hessians) | convex problems | yes | no | start point |
| nlpbb | nloptr | MINLP, MIQCQP, MIQP | convex relaxations | no | no | no |
| roi | ROI + plugin | LP, MILP, convex QP (linear rows only) | LP, MILP | no | no | no |

HiGHS is required; the others are optional. The IPOPT backend is compiled
in when `pkg-config` finds IPOPT at installation (`brew install ipopt`,
`apt install coinor-libipopt-dev`); check with `opt_ipopt_available()`.

## Installation

```r
# install.packages("remotes")
remotes::install_github("adamrsophie/optcompR")
```

Optional integrations: `remotes::install_github("adamrsophie/z3r")` (needs
the Z3 library), `remotes::install_github("adamrsophie/symcompR")`,
`remotes::install_github("adamrsophie/jaxR")`,
`remotes::install_github("adamrsophie/modelcompR")`.

## Documentation

Articles at <https://adamrsophie.github.io/optcompR/> cover linear,
quadratic, mixed-integer, nonlinear, conic, robust, stochastic and
multiobjective optimization, indicator and SOS constraints, parameterized and
sparse models, solver comparison, diagnostics, infeasible models,
sensitivity, derivatives of optimal solutions, exact certificates, ODE
fitting and the underlying concepts.

## License

GPL-3.
