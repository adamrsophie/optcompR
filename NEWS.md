# optcompR 0.1.0

First release. optcompR is a solver independent modeling layer for
mathematical optimization in R: models are written with ordinary R
operators, classified with recorded evidence, lowered to sparse standard
forms, solved by installed backends, and returned as standardized,
independently checked results.

## Modeling

* Continuous, integer and binary decision variables, parameters and
  symmetric matrix variables (`opt_variable()`, `opt_parameter()`,
  `opt_sym_variable()`), with vectorized linear, quadratic and nonlinear
  expressions, sparse `Matrix` support, ranged constraints and named
  blocks. Matrix products follow the conformability rules of base R.
* Cone constraints: second-order, exponential, power and positive
  semidefinite (`opt_soc_constraint()`, `opt_exp_cone()`, `opt_pow_cone()`,
  `opt_psd_constraint()`).
* Indicator, SOS1, SOS2, complementarity and logical constraints
  (`opt_indicator()`, `opt_sos1()`, `opt_sos2()`, `opt_implies()`,
  `opt_and()`, `opt_or()`), reformulated exactly for solvers without native
  support.
* External R functions with supplied, jaxR or labeled finite difference
  gradients (`opt_function()`).
* Model inspection and editing: `opt_copy()`, `opt_variables()`,
  `opt_constraints()`, `opt_special_constraints()`,
  `opt_delete_constraint()`, `opt_terms()`, `opt_degree()`, `opt_symbols()`,
  `opt_evaluate()`, `opt_standard_form()` and `opt_row_hessian()`.

## Classification and solver selection

* `opt_classify()` assigns the class (LP, MILP, QP, MIQP, QCQP, MIQCQP,
  SOCP, MISOCP, CONIC, NLP, MINLP) and records how convexity was
  established: structurally, numerically, or by disciplined convex
  programming rules with interval bounds for nonlinear expressions. When
  convexity cannot be proven the model is reported as `"unknown"`, never
  guessed.
* `opt_solve()` selects an installed backend whose R binding supports the
  class, preferring backends with global guarantees for nonconvex models. A
  model is never approximated to fit a solver; unsupported models are
  refused with the reason.

## Backends

* HiGHS (LP, MILP, convex QP, sensitivity ranging, warm starts and in-place
  updates across repeated solves), SCIP (global MIQCQP, MISOCP and
  polynomial NLP and MINLP by exact lifting to quadratic rows, native
  indicator and SOS constraints), Clarabel and SCS (LP, convex QP and conic
  programs, with cone multipliers), OSQP (convex QP), NLopt (local NLP),
  IPOPT with exact sparse Hessians (compiled in when available), NLP branch
  and bound for MINLP with outer approximation bounds certified by weak
  duality (global when the continuous relaxation is convex, a heuristic
  otherwise), and ROI plugins. A local solver that stops because its steps
  stalled reports a local optimum, even on a convex model.
* `opt_register_solver()` makes the backend contract public, so other
  solvers can be connected without changing optcompR.

## Results and diagnostics

* Standardized termination and point statuses that separate proven
  optimality, local optimality, early stops (including `"gap_limit"` with
  its gap and a valid bound) and infeasibility certificates.
* Every returned point is checked against the model before it is reported
  feasible, at tolerances that follow the accuracy requested from the
  solver and are stored in the result; a point where the objective is not
  defined is not a solution.
* Duals are derivatives of the optimal value with a documented sign
  convention; `opt_check()` recomputes feasibility and the Karush Kuhn
  Tucker conditions, including conic ones, from the original expressions.
* Exact gradients, Jacobians and Lagrangian Hessians.
* `opt_explain_infeasibility()` reports which rows and bounds must move,
  and by how much, through a weighted elastic (L1) relaxation.
* Exact feasibility checks and optimality certificates for LP and MILP
  through z3r (`opt_z3_check()`, `opt_z3_certify()`).

## Sensitivity and analysis

* Sensitivity ranging for LPs (`opt_sensitivity()`), reporting the active
  side of ranged rows and covering problems with bounds only, fixed relaxations of
  mixed-integer results (`opt_fixed_relaxation()`) and parameter scenario
  sweeps (`opt_scenarios()`).
* Derivatives of optimal solutions with respect to parameters
  (`opt_solution_jacobian()`), for smooth, conic and semidefinite programs.
* Multiobjective scalarizations: `opt_weighted()`, `opt_lexicographic()`
  and `opt_epsilon_constraint()`.

## Robust and stochastic optimization

* `opt_robust_constraint()` (box, ellipsoidal and budget uncertainty),
  `opt_chance_constraint()`, `opt_cvar()`, `opt_expectation()` and
  `opt_worst_expectation()`, all as exact reformulations.

## Integrations

* ODE least squares losses for modelcompR models with forward or adjoint
  gradients (`opt_ode_loss()`), both covering parameter-dependent initial
  values, events and observations in any order.
* Conversion to and from symcompR expressions, to CVXR problems
  (`opt_to_cvxr()`), and export to MPS and LP files.

## Documentation

* Twenty vignettes and a pkgdown site.
