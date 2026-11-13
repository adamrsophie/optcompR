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
  IPOPT with exact Hessians (compiled in when found with pkg-config),
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
* Correctness safeguards: binary variables keep 0/1 bounds whatever
  `opt_set_bounds()` sets; a smallest Hessian eigenvalue that is negative
  beyond rounding but inside the curvature tolerance makes convexity
  "unknown", and convex-only solvers now require proven convexity; NLopt's
  objective-target stop is `"objective_limit"`; OSQP and IPOPT points at a
  limit are reported feasible only after a feasibility check; IPOPT
  diverging iterates are no longer called unbounded; `nlpbb` no longer fails
  on its first node for convex models; SOS2 checks catch two nonadjacent
  nonzeros; constraints with a lower bound of `Inf` or an upper bound of
  `-Inf` are rejected.
* `log(x, base)` honors `base`; powers with exponents 0 and 1 have exact
  derivatives at zero; `t(x) %*% (A %*% x)` works.
* `opt_solution_jacobian()` requires an optimal termination, the KKT
  conditions and the second-order sufficient condition, and handles
  solutions where every variable sits at a bound. Estimated multipliers
  respect sign restrictions.
* `opt_fixed_relaxation()` resolves indicator and SOS constraints at the
  incumbent. Cone and indicator constructors accept products with
  parameters. Constraint references from another model are refused.
  `opt_scenarios()` validates values. `warm_start = FALSE` is honored by
  NLopt and IPOPT. IPOPT options are type checked and rejected options are
  errors. Exact Z3 coefficients cover subnormal numbers.
* Local NLP solvers start strictly inside the bounds when no start is
  given, and stop with a clear message when a function is not finite at the
  start. NLopt points with NaN are reported as no solution.
* `nlpbb` counts nodes without a certified verdict as lost, reports
  `"unknown"` instead of infeasible when nodes were lost, catches node
  errors, and re-solves with integers fixed before accepting a rounded
  incumbent.
* Clarabel accepts quadratic objectives with cross terms and returns cone
  multipliers for models with cones only. ROI loads its plugins itself and
  passes equality rows correctly.
* IPOPT: single precision builds fall back to the stub, older IPOPT uses
  `max_cpu_time` for time limits, options are tried in the type of the R
  value first, and callback errors keep their message; a failed trial
  evaluation is a warning when IPOPT recovers.
* Curvature tests are relative to the largest eigenvalue and the large
  matrix Cholesky test shifts only at rounding level. Quadratic rows
  without finite bounds no longer make a model nonconvex.
* Matrix indices `X[cbind(i, j)]` select coordinates. `mean()` works on
  expressions. Element names must be unique across blocks.
* `opt_solution_jacobian()` checks primal feasibility, refuses `abs()` at
  a kink, and treats a point near a bound with a zero multiplier as
  interior. Row and bound feasibility scale each side by its own bound.
* Indicator switches must stay binary. `opt_fixed_relaxation()` refuses a
  result whose model changed since the solve. Automatic selection skips
  solvers that cannot reformulate indicator or SOS constraints, and
  controls a backend ignores are reported in the result warnings.
* `opt_z3_certify()` also proves that the reported value is attained, and
  `opt_z3_check()` works with nonlinear objectives.
