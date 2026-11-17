# Changelog

## optcompR 0.1.0

First release. optcompR is a solver independent modeling layer for
mathematical optimization in R: models are written with ordinary R
operators, classified with recorded evidence, lowered to sparse standard
forms, solved by installed backends, and returned as standardized,
independently checked results.

### Modeling

- Continuous, integer and binary decision variables, parameters and
  symmetric matrix variables
  ([`opt_variable()`](https://adamrsophie.github.io/optcompR/reference/opt_variable.md),
  [`opt_parameter()`](https://adamrsophie.github.io/optcompR/reference/opt_parameter.md),
  [`opt_sym_variable()`](https://adamrsophie.github.io/optcompR/reference/opt_sym_variable.md)),
  with vectorized linear, quadratic and nonlinear expressions, sparse
  `Matrix` support, ranged constraints and named blocks. Matrix products
  follow the conformability rules of base R.
- Cone constraints: second-order, exponential, power and positive
  semidefinite
  ([`opt_soc_constraint()`](https://adamrsophie.github.io/optcompR/reference/opt_soc_constraint.md),
  [`opt_exp_cone()`](https://adamrsophie.github.io/optcompR/reference/opt_soc_constraint.md),
  [`opt_pow_cone()`](https://adamrsophie.github.io/optcompR/reference/opt_soc_constraint.md),
  [`opt_psd_constraint()`](https://adamrsophie.github.io/optcompR/reference/opt_soc_constraint.md)).
- Indicator, SOS1, SOS2, complementarity and logical constraints
  ([`opt_indicator()`](https://adamrsophie.github.io/optcompR/reference/opt_indicator.md),
  [`opt_sos1()`](https://adamrsophie.github.io/optcompR/reference/opt_sos1.md),
  [`opt_sos2()`](https://adamrsophie.github.io/optcompR/reference/opt_sos1.md),
  [`opt_implies()`](https://adamrsophie.github.io/optcompR/reference/opt_implies.md),
  [`opt_and()`](https://adamrsophie.github.io/optcompR/reference/opt_implies.md),
  [`opt_or()`](https://adamrsophie.github.io/optcompR/reference/opt_implies.md)),
  reformulated exactly for solvers without native support.
- External R functions with supplied, jaxR or labeled finite difference
  gradients
  ([`opt_function()`](https://adamrsophie.github.io/optcompR/reference/opt_function.md)).
- Model inspection and editing:
  [`opt_copy()`](https://adamrsophie.github.io/optcompR/reference/opt_copy.md),
  [`opt_variables()`](https://adamrsophie.github.io/optcompR/reference/opt_variables.md),
  [`opt_constraints()`](https://adamrsophie.github.io/optcompR/reference/opt_constraints.md),
  [`opt_special_constraints()`](https://adamrsophie.github.io/optcompR/reference/opt_special_constraints.md),
  [`opt_delete_constraint()`](https://adamrsophie.github.io/optcompR/reference/opt_delete_constraint.md),
  [`opt_terms()`](https://adamrsophie.github.io/optcompR/reference/opt_terms.md),
  [`opt_degree()`](https://adamrsophie.github.io/optcompR/reference/opt_terms.md),
  [`opt_symbols()`](https://adamrsophie.github.io/optcompR/reference/opt_terms.md),
  [`opt_evaluate()`](https://adamrsophie.github.io/optcompR/reference/opt_evaluate.md),
  [`opt_standard_form()`](https://adamrsophie.github.io/optcompR/reference/opt_standard_form.md)
  and
  [`opt_row_hessian()`](https://adamrsophie.github.io/optcompR/reference/opt_row_hessian.md).

### Classification and solver selection

- [`opt_classify()`](https://adamrsophie.github.io/optcompR/reference/opt_classify.md)
  assigns the class (LP, MILP, QP, MIQP, QCQP, MIQCQP, SOCP, MISOCP,
  CONIC, NLP, MINLP) and records how convexity was established:
  structurally, numerically, or by disciplined convex programming rules
  with interval bounds for nonlinear expressions. When convexity cannot
  be proven the model is reported as `"unknown"`, never guessed.
- [`opt_solve()`](https://adamrsophie.github.io/optcompR/reference/opt_solve.md)
  selects an installed backend whose R binding supports the class,
  preferring backends with global guarantees for nonconvex models. A
  model is never approximated to fit a solver; unsupported models are
  refused with the reason.

### Backends

- HiGHS (LP, MILP, convex QP, sensitivity ranging, warm starts and
  in-place updates across repeated solves), SCIP (global MIQCQP, MISOCP
  and polynomial NLP and MINLP by exact lifting to quadratic rows,
  native indicator and SOS constraints), Clarabel and SCS (LP, convex QP
  and conic programs, with cone multipliers), OSQP (convex QP), NLopt
  (local NLP), IPOPT with exact sparse Hessians (compiled in when
  available), NLP branch and bound for MINLP with outer approximation
  bounds certified by weak duality, and ROI plugins.
- [`opt_register_solver()`](https://adamrsophie.github.io/optcompR/reference/opt_register_solver.md)
  makes the backend contract public, so other solvers can be connected
  without changing optcompR.

### Results and diagnostics

- Standardized termination and point statuses that separate proven
  optimality, local optimality, early stops (including `"gap_limit"`
  with its gap and a valid bound) and infeasibility certificates.
- Every returned point is checked against the model before it is
  reported feasible, at tolerances that follow the accuracy requested
  from the solver and are stored in the result; a point where the
  objective is not defined is not a solution.
- Duals are derivatives of the optimal value with a documented sign
  convention;
  [`opt_check()`](https://adamrsophie.github.io/optcompR/reference/opt_check.md)
  recomputes feasibility and the Karush Kuhn Tucker conditions,
  including conic ones, from the original expressions.
- Exact gradients, Jacobians and Lagrangian Hessians.
- [`opt_explain_infeasibility()`](https://adamrsophie.github.io/optcompR/reference/opt_explain_infeasibility.md)
  reports which rows and bounds must move, and by how much, through a
  weighted elastic (L1) relaxation.
- Exact feasibility checks and optimality certificates for LP and MILP
  through z3r
  ([`opt_z3_check()`](https://adamrsophie.github.io/optcompR/reference/opt_z3_check.md),
  [`opt_z3_certify()`](https://adamrsophie.github.io/optcompR/reference/opt_z3_certify.md)).

### Sensitivity and analysis

- Sensitivity ranging for LPs
  ([`opt_sensitivity()`](https://adamrsophie.github.io/optcompR/reference/opt_sensitivity.md)),
  fixed relaxations of mixed-integer results
  ([`opt_fixed_relaxation()`](https://adamrsophie.github.io/optcompR/reference/opt_fixed_relaxation.md))
  and parameter scenario sweeps
  ([`opt_scenarios()`](https://adamrsophie.github.io/optcompR/reference/opt_scenarios.md)).
- Derivatives of optimal solutions with respect to parameters
  ([`opt_solution_jacobian()`](https://adamrsophie.github.io/optcompR/reference/opt_solution_jacobian.md)),
  for smooth, conic and semidefinite programs.
- Multiobjective scalarizations:
  [`opt_weighted()`](https://adamrsophie.github.io/optcompR/reference/opt_weighted.md),
  [`opt_lexicographic()`](https://adamrsophie.github.io/optcompR/reference/opt_lexicographic.md)
  and
  [`opt_epsilon_constraint()`](https://adamrsophie.github.io/optcompR/reference/opt_epsilon_constraint.md).

### Robust and stochastic optimization

- [`opt_robust_constraint()`](https://adamrsophie.github.io/optcompR/reference/opt_robust_constraint.md)
  (box, ellipsoidal and budget uncertainty),
  [`opt_chance_constraint()`](https://adamrsophie.github.io/optcompR/reference/opt_chance_constraint.md),
  [`opt_cvar()`](https://adamrsophie.github.io/optcompR/reference/opt_cvar.md),
  [`opt_expectation()`](https://adamrsophie.github.io/optcompR/reference/opt_cvar.md)
  and
  [`opt_worst_expectation()`](https://adamrsophie.github.io/optcompR/reference/opt_cvar.md),
  all as exact reformulations.

### Integrations

- ODE least squares losses for modelcompR models with forward or adjoint
  gradients
  ([`opt_ode_loss()`](https://adamrsophie.github.io/optcompR/reference/opt_ode_loss.md)),
  both covering parameter-dependent initial values, events and
  observations in any order.
- Conversion to and from symcompR expressions, to CVXR problems
  ([`opt_to_cvxr()`](https://adamrsophie.github.io/optcompR/reference/opt_to_cvxr.md)),
  and export to MPS and LP files.

### Documentation

- Twenty vignettes and a pkgdown site.
