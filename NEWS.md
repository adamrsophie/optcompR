# optcompR 0.1.0

* First release.
* Modeling: continuous, integer and binary variables, parameters, vectorized
  linear, quadratic and nonlinear expressions with sparse `Matrix` support,
  ranged constraints and named blocks.
* Classification into LP, MILP, QP, MIQP, QCQP, MIQCQP, SOCP, MISOCP,
  CONIC, NLP and MINLP with recorded convexity evidence, and
  capability-aware solver selection.
* Backends: HiGHS (LP, MILP, convex QP, ranging, warm starts), SCIP (global
  MIQCQP and MISOCP, native indicator and SOS constraints), Clarabel (LP,
  convex QP and conic programs), OSQP (convex QP), NLopt SLSQP (local NLP),
  NLP branch and bound for MINLP, and ROI plugins.
* Cone constraints (`opt_soc_constraint()`, `opt_exp_cone()`,
  `opt_pow_cone()`, `opt_psd_constraint()`) and symmetric matrix variables
  (`opt_sym_variable()`).
* Indicator, SOS1, SOS2, complementarity and logical constraints with exact
  reformulations for solvers without native support.
* Robust and stochastic optimization: `opt_robust_constraint()`,
  `opt_chance_constraint()`, `opt_cvar()`, `opt_expectation()`,
  `opt_worst_expectation()`.
* Derivatives of optimal solutions with respect to parameters
  (`opt_solution_jacobian()`).
* ODE least squares losses for modelcompR models (`opt_ode_loss()`).
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
* Twenty vignettes and a pkgdown site.
