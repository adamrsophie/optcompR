# Changelog

## optcompR 0.1.0

- First release.
- Modeling: continuous, integer and binary variables, parameters,
  vectorized linear, quadratic and nonlinear expressions with sparse
  `Matrix` support, ranged constraints and named blocks.
- Classification into LP, MILP, QP, MIQP, QCQP, MIQCQP, NLP and MINLP
  with recorded convexity evidence, and capability-aware solver
  selection.
- Backends: HiGHS (LP, MILP, convex QP, ranging, warm starts), SCIP
  (global MIQCQP), OSQP (convex QP), NLopt SLSQP (local NLP) and ROI
  plugins.
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
- Fifteen vignettes and a pkgdown site.
