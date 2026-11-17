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
  SLSQP (local NLP), IPOPT with exact Hessians (compiled in when found
  with pkg-config), NLP branch and bound for MINLP, and ROI plugins.
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
- Correctness safeguards: binary variables keep 0/1 bounds whatever
  [`opt_set_bounds()`](https://adamrsophie.github.io/optcompR/reference/opt_set_bounds.md)
  sets; a smallest Hessian eigenvalue that is negative beyond rounding
  but inside the curvature tolerance makes convexity “unknown”, and
  convex-only solvers now require proven convexity; NLopt’s
  objective-target stop is `"objective_limit"`; OSQP and IPOPT points at
  a limit are reported feasible only after a feasibility check; IPOPT
  diverging iterates are no longer called unbounded; `nlpbb` no longer
  fails on its first node for convex models; SOS2 checks catch two
  nonadjacent nonzeros; constraints with a lower bound of `Inf` or an
  upper bound of `-Inf` are rejected.
- `log(x, base)` honors `base`; powers with exponents 0 and 1 have exact
  derivatives at zero; `t(x) %*% (A %*% x)` works.
- [`opt_solution_jacobian()`](https://adamrsophie.github.io/optcompR/reference/opt_solution_jacobian.md)
  requires an optimal termination, the KKT conditions and the
  second-order sufficient condition, and handles solutions where every
  variable sits at a bound. Estimated multipliers respect sign
  restrictions.
- [`opt_fixed_relaxation()`](https://adamrsophie.github.io/optcompR/reference/opt_fixed_relaxation.md)
  resolves indicator and SOS constraints at the incumbent. Cone and
  indicator constructors accept products with parameters. Constraint
  references from another model are refused.
  [`opt_scenarios()`](https://adamrsophie.github.io/optcompR/reference/opt_scenarios.md)
  validates values. `warm_start = FALSE` is honored by NLopt and IPOPT.
  IPOPT options are type checked and rejected options are errors. Exact
  Z3 coefficients cover subnormal numbers.
- Local NLP solvers start strictly inside the bounds when no start is
  given, and stop with a clear message when a function is not finite at
  the start. NLopt points with NaN are reported as no solution.
- `nlpbb` counts nodes without a certified verdict as lost, reports
  `"unknown"` instead of infeasible when nodes were lost, catches node
  errors, and re-solves with integers fixed before accepting a rounded
  incumbent.
- Clarabel accepts quadratic objectives with cross terms and returns
  cone multipliers for models with cones only. ROI loads its plugins
  itself and passes equality rows correctly.
- IPOPT: single precision builds fall back to the stub, older IPOPT uses
  `max_cpu_time` for time limits, options are tried in the type of the R
  value first, and callback errors keep their message; a failed trial
  evaluation is a warning when IPOPT recovers.
- Curvature tests are relative to the largest eigenvalue and the large
  matrix Cholesky test shifts only at rounding level. Quadratic rows
  without finite bounds no longer make a model nonconvex.
- Matrix indices `X[cbind(i, j)]` select coordinates.
  [`mean()`](https://rdrr.io/r/base/mean.html) works on expressions.
  Element names must be unique across blocks.
- [`opt_solution_jacobian()`](https://adamrsophie.github.io/optcompR/reference/opt_solution_jacobian.md)
  checks primal feasibility, refuses
  [`abs()`](https://rdrr.io/r/base/MathFun.html) at a kink, and treats a
  point near a bound with a zero multiplier as interior. Row and bound
  feasibility scale each side by its own bound.
- Indicator switches must stay binary.
  [`opt_fixed_relaxation()`](https://adamrsophie.github.io/optcompR/reference/opt_fixed_relaxation.md)
  refuses a result whose model changed since the solve. Automatic
  selection skips solvers that cannot reformulate indicator or SOS
  constraints, and controls a backend ignores are reported in the result
  warnings.
- [`opt_z3_certify()`](https://adamrsophie.github.io/optcompR/reference/opt_z3_certify.md)
  also proves that the reported value is attained, and
  [`opt_z3_check()`](https://adamrsophie.github.io/optcompR/reference/opt_z3_check.md)
  works with nonlinear objectives.
- A power whose exponent is a parameter expression stays a power node
  and becomes a constant power once values are substituted, so `x^p`
  with `p = 2` is a quadratic even for negative `x`; derivatives in the
  exponent are exact for positive bases and refused where they do not
  exist.
- [`opt_copy()`](https://adamrsophie.github.io/optcompR/reference/opt_copy.md)
  copies a model.
  [`opt_special_constraints()`](https://adamrsophie.github.io/optcompR/reference/opt_special_constraints.md)
  lists cones, indicator constraints and special ordered sets, and
  [`opt_delete_constraint()`](https://adamrsophie.github.io/optcompR/reference/opt_delete_constraint.md)
  deletes them.
  [`opt_terms()`](https://adamrsophie.github.io/optcompR/reference/opt_terms.md),
  [`opt_degree()`](https://adamrsophie.github.io/optcompR/reference/opt_terms.md),
  [`opt_symbols()`](https://adamrsophie.github.io/optcompR/reference/opt_terms.md)
  and
  [`opt_evaluate()`](https://adamrsophie.github.io/optcompR/reference/opt_evaluate.md)
  expose the structure and values of expressions.
- A mixed-integer solve that stops at a user `mip_gap` with a positive
  gap reports `"gap_limit"` instead of `"optimal"`. SCIP results carry
  the dual bound recovered from SCIP’s gap.
- `configure` links a test program before enabling IPOPT, and Windows
  builds use IPOPT from a COIN-OR binary release named by `IPOPT_HOME`.
- [`opt_register_solver()`](https://adamrsophie.github.io/optcompR/reference/opt_register_solver.md)
  and
  [`opt_unregister_solver()`](https://adamrsophie.github.io/optcompR/reference/opt_register_solver.md)
  make the backend contract public; what a backend returns is validated.
- SCS is a second conic backend (`solver = "scs"`). Clarabel and SCS
  results carry cone multipliers by cone in `extra$cone_dual`, and
  [`opt_check()`](https://adamrsophie.github.io/optcompR/reference/opt_check.md)
  verifies the conic KKT conditions with them (stationarity, dual cone
  membership and complementarity per cone).
- [`opt_solution_jacobian()`](https://adamrsophie.github.io/optcompR/reference/opt_solution_jacobian.md)
  differentiates conic solutions through the smooth boundaries of active
  second-order, exponential and power cones, and through the linearized
  complementarity of active PSD constraints; inactive cones drop out.
- Lowered constraint blocks and objectives are reused across parameter
  updates when the parameters they use keep their values, and HiGHS
  restarts LPs from the previous optimal basis. Scalar constraints built
  in loops are about a quarter faster.
- SCIP solves polynomial NLP and MINLP models to global optimality:
  products of degree above two are lifted exactly to quadratic rows with
  auxiliary variables bounded by interval arithmetic. Nonconvex
  polynomial models choose SCIP automatically; models with other
  operators are refused with the reason.
- [`opt_to_cvxr()`](https://adamrsophie.github.io/optcompR/reference/opt_to_cvxr.md)
  converts linear, quadratic, mixed-integer and conic models to CVXR
  problems, and
  [`opt_cvxr_values()`](https://adamrsophie.github.io/optcompR/reference/opt_cvxr_values.md)
  reads the solution back by variable name.
- [`opt_explain_infeasibility()`](https://adamrsophie.github.io/optcompR/reference/opt_explain_infeasibility.md)
  solves a weighted elastic relaxation of the rows and variable bounds
  and reports which must move, and by how much, for an infeasible model
  to become feasible.
- [`opt_classify()`](https://adamrsophie.github.io/optcompR/reference/opt_classify.md)
  proves convexity of nonlinear models by composition rules (disciplined
  convex programming) with interval bounds from the variable bounds, so
  local solvers and NLP branch and bound report global optima for more
  models.
- The IPOPT binding passes sparse Jacobian and Hessian structures
  derived from the expressions instead of dense ones.
