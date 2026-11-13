# Changelog

## optcompR 0.1.0

- First release.
- Modeling: continuous, integer and binary variables, parameters,
  vectorized linear, quadratic and nonlinear expressions with sparse
  `Matrix` support, ranged constraints and named blocks.
- Classification into LP, MILP, QP, MIQP, QCQP, MIQCQP, SOCP, MISOCP,
  CONIC, NLP and MINLP with recorded convexity evidence, and
  capability-aware solver selection.
- Backends: HiGHS (LP, MILP, convex QP, ranging, warm starts), SCIP
  (global MIQCQP and MISOCP, native indicator and SOS constraints),
  Clarabel (LP, convex QP and conic programs), OSQP (convex QP), NLopt
  SLSQP (local NLP), NLP branch and bound for MINLP, and ROI plugins.
- Cone constraints
  ([`opt_soc_constraint()`](https://adamrsophie.github.io/optcompR/reference/opt_soc_constraint.md),
  [`opt_exp_cone()`](https://adamrsophie.github.io/optcompR/reference/opt_soc_constraint.md),
  [`opt_pow_cone()`](https://adamrsophie.github.io/optcompR/reference/opt_soc_constraint.md),
  [`opt_psd_constraint()`](https://adamrsophie.github.io/optcompR/reference/opt_soc_constraint.md))
  and symmetric matrix variables
  ([`opt_sym_variable()`](https://adamrsophie.github.io/optcompR/reference/opt_sym_variable.md)).
- Indicator, SOS1, SOS2, complementarity and logical constraints with
  exact reformulations for solvers without native support.
- Robust and stochastic optimization:
  [`opt_robust_constraint()`](https://adamrsophie.github.io/optcompR/reference/opt_robust_constraint.md),
  [`opt_chance_constraint()`](https://adamrsophie.github.io/optcompR/reference/opt_chance_constraint.md),
  [`opt_cvar()`](https://adamrsophie.github.io/optcompR/reference/opt_cvar.md),
  [`opt_expectation()`](https://adamrsophie.github.io/optcompR/reference/opt_cvar.md),
  [`opt_worst_expectation()`](https://adamrsophie.github.io/optcompR/reference/opt_cvar.md).
- Derivatives of optimal solutions with respect to parameters
  ([`opt_solution_jacobian()`](https://adamrsophie.github.io/optcompR/reference/opt_solution_jacobian.md)).
- ODE least squares losses for modelcompR models
  ([`opt_ode_loss()`](https://adamrsophie.github.io/optcompR/reference/opt_ode_loss.md)).
- Results with standardized statuses, duals as derivatives of the
  optimal value, independent feasibility and KKT diagnostics
  ([`opt_check()`](https://adamrsophie.github.io/optcompR/reference/opt_check.md)).
- Exact gradients, Jacobians and Lagrangian Hessians; external R
  functions with supplied, jaxR or labeled finite difference gradients
  ([`opt_function()`](https://adamrsophie.github.io/optcompR/reference/opt_function.md)).
- Sensitivity ranging, fixed relaxations of MILPs and scenario sweeps.
- Multiobjective scalarizations:
  [`opt_weighted()`](https://adamrsophie.github.io/optcompR/reference/opt_weighted.md),
  [`opt_lexicographic()`](https://adamrsophie.github.io/optcompR/reference/opt_lexicographic.md),
  [`opt_epsilon_constraint()`](https://adamrsophie.github.io/optcompR/reference/opt_epsilon_constraint.md).
- Exact feasibility checks and optimality certificates for LP and MILP
  through z3r
  ([`opt_z3_check()`](https://adamrsophie.github.io/optcompR/reference/opt_z3_check.md),
  [`opt_z3_certify()`](https://adamrsophie.github.io/optcompR/reference/opt_z3_certify.md)).
- Conversion to and from symcompR expressions; MPS and LP export.
- Twenty vignettes and a pkgdown site.
