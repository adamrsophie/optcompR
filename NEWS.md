# optcompR 0.1.0

* First release.
* Modeling: continuous, integer and binary variables, parameters, vectorized
  linear, quadratic and nonlinear expressions with sparse `Matrix` support,
  ranged constraints and named blocks.
* Classification into LP, MILP, QP, MIQP, QCQP, MIQCQP, NLP and MINLP with
  recorded convexity evidence, and capability-aware solver selection.
* Backends: HiGHS (LP, MILP, convex QP, ranging, warm starts), SCIP (global
  MIQCQP), OSQP (convex QP), NLopt SLSQP (local NLP) and ROI plugins.
* Results with standardized statuses, duals as derivatives of the optimal
  value, independent feasibility and KKT diagnostics (`opt_check()`).
* Exact gradients, Jacobians and Lagrangian Hessians; external R functions
  with supplied, jaxR or labeled finite difference gradients
  (`opt_function()`).
* Sensitivity ranging, fixed relaxations of MILPs and scenario sweeps.
* Multiobjective scalarizations: `opt_weighted()`, `opt_lexicographic()`,
  `opt_epsilon_constraint()`.
* Exact feasibility checks and optimality certificates for LP and MILP
  through z3r (`opt_z3_check()`, `opt_z3_certify()`).
* Conversion to and from symcompR expressions; MPS and LP export.
* Fifteen vignettes and a pkgdown site.
