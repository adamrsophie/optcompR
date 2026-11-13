# Package index

## Models and variables

- [`opt_model()`](https://adamrsophie.github.io/optcompR/reference/opt_model.md)
  : Create an optimization model
- [`opt_variable()`](https://adamrsophie.github.io/optcompR/reference/opt_variable.md)
  : Add decision variables to a model
- [`opt_parameter()`](https://adamrsophie.github.io/optcompR/reference/opt_parameter.md)
  : Add parameters to a model
- [`opt_get()`](https://adamrsophie.github.io/optcompR/reference/opt_get.md)
  : Retrieve a variable or parameter block by name
- [`opt_variables()`](https://adamrsophie.github.io/optcompR/reference/opt_variables.md)
  : List the variables of a model
- [`opt_set_bounds()`](https://adamrsophie.github.io/optcompR/reference/opt_set_bounds.md)
  : Change variable bounds
- [`opt_set_type()`](https://adamrsophie.github.io/optcompR/reference/opt_set_type.md)
  : Change variable types
- [`opt_fix()`](https://adamrsophie.github.io/optcompR/reference/opt_fix.md)
  : Fix variables to values
- [`opt_set_start()`](https://adamrsophie.github.io/optcompR/reference/opt_set_start.md)
  : Set start values for variables
- [`opt_set_parameter()`](https://adamrsophie.github.io/optcompR/reference/opt_set_parameter.md)
  : Change parameter values

## Expressions

- [`Ops(`*`<opt_expr>`*`)`](https://adamrsophie.github.io/optcompR/reference/opt_operators.md)
  [`Math(`*`<opt_expr>`*`)`](https://adamrsophie.github.io/optcompR/reference/opt_operators.md)
  [`Summary(`*`<opt_expr>`*`)`](https://adamrsophie.github.io/optcompR/reference/opt_operators.md)
  [`mean(`*`<opt_expr>`*`)`](https://adamrsophie.github.io/optcompR/reference/opt_operators.md)
  [`matrixOps(`*`<opt_expr>`*`)`](https://adamrsophie.github.io/optcompR/reference/opt_operators.md)
  [`` `%*%`( ``*`<Matrix>`*`,`*`<opt_expr>`*`)`](https://adamrsophie.github.io/optcompR/reference/opt_operators.md)
  [`` `%*%`( ``*`<opt_expr>`*`,`*`<Matrix>`*`)`](https://adamrsophie.github.io/optcompR/reference/opt_operators.md)
  : Arithmetic and comparison operators for expressions
- [`` `[`( ``*`<opt_expr>`*`)`](https://adamrsophie.github.io/optcompR/reference/opt_expr_shape.md)
  [`` `[[`( ``*`<opt_expr>`*`)`](https://adamrsophie.github.io/optcompR/reference/opt_expr_shape.md)
  [`length(`*`<opt_expr>`*`)`](https://adamrsophie.github.io/optcompR/reference/opt_expr_shape.md)
  [`dim(`*`<opt_expr>`*`)`](https://adamrsophie.github.io/optcompR/reference/opt_expr_shape.md)
  [`` `dim<-`( ``*`<opt_expr>`*`)`](https://adamrsophie.github.io/optcompR/reference/opt_expr_shape.md)
  [`t(`*`<opt_expr>`*`)`](https://adamrsophie.github.io/optcompR/reference/opt_expr_shape.md)
  [`c(`*`<opt_expr>`*`)`](https://adamrsophie.github.io/optcompR/reference/opt_expr_shape.md)
  : Shape and indexing of expressions
- [`is_opt_expr()`](https://adamrsophie.github.io/optcompR/reference/is_opt_expr.md)
  : Test whether an object is an optimization expression
- [`opt_function()`](https://adamrsophie.github.io/optcompR/reference/opt_function.md)
  : Wrap an R function as a nonlinear expression

## Constraints and objectives

- [`opt_constraint()`](https://adamrsophie.github.io/optcompR/reference/opt_constraint.md)
  : Add constraints to a model
- [`opt_constraints()`](https://adamrsophie.github.io/optcompR/reference/opt_constraints.md)
  : List the constraint blocks of a model
- [`opt_delete_constraint()`](https://adamrsophie.github.io/optcompR/reference/opt_delete_constraint.md)
  : Delete a constraint block
- [`opt_objective()`](https://adamrsophie.github.io/optcompR/reference/opt_objective.md)
  : Set the objective of a model

## Cones, discrete structure and logic

- [`opt_soc_constraint()`](https://adamrsophie.github.io/optcompR/reference/opt_soc_constraint.md)
  [`opt_exp_cone()`](https://adamrsophie.github.io/optcompR/reference/opt_soc_constraint.md)
  [`opt_pow_cone()`](https://adamrsophie.github.io/optcompR/reference/opt_soc_constraint.md)
  [`opt_psd_constraint()`](https://adamrsophie.github.io/optcompR/reference/opt_soc_constraint.md)
  : Conic constraints
- [`opt_sym_variable()`](https://adamrsophie.github.io/optcompR/reference/opt_sym_variable.md)
  : Symmetric matrix variables
- [`opt_indicator()`](https://adamrsophie.github.io/optcompR/reference/opt_indicator.md)
  : Indicator constraints
- [`opt_sos1()`](https://adamrsophie.github.io/optcompR/reference/opt_sos1.md)
  [`opt_sos2()`](https://adamrsophie.github.io/optcompR/reference/opt_sos1.md)
  [`opt_complementarity()`](https://adamrsophie.github.io/optcompR/reference/opt_sos1.md)
  : Special ordered sets and complementarity
- [`opt_implies()`](https://adamrsophie.github.io/optcompR/reference/opt_implies.md)
  [`opt_and()`](https://adamrsophie.github.io/optcompR/reference/opt_implies.md)
  [`opt_or()`](https://adamrsophie.github.io/optcompR/reference/opt_implies.md)
  : Logical constraints on binary variables

## Robust and stochastic optimization

- [`opt_robust_constraint()`](https://adamrsophie.github.io/optcompR/reference/opt_robust_constraint.md)
  : Robust linear constraints
- [`opt_chance_constraint()`](https://adamrsophie.github.io/optcompR/reference/opt_chance_constraint.md)
  : Chance constraints
- [`opt_cvar()`](https://adamrsophie.github.io/optcompR/reference/opt_cvar.md)
  [`opt_expectation()`](https://adamrsophie.github.io/optcompR/reference/opt_cvar.md)
  [`opt_worst_expectation()`](https://adamrsophie.github.io/optcompR/reference/opt_cvar.md)
  : Risk measures and expectations over scenarios

## Solving

- [`opt_solve()`](https://adamrsophie.github.io/optcompR/reference/opt_solve.md)
  : Solve an optimization model
- [`opt_solvers()`](https://adamrsophie.github.io/optcompR/reference/opt_solvers.md)
  : Available optimization solvers
- [`opt_solver_capabilities()`](https://adamrsophie.github.io/optcompR/reference/opt_solver_capabilities.md)
  : Capabilities of one solver backend
- [`opt_ipopt_available()`](https://adamrsophie.github.io/optcompR/reference/opt_ipopt_available.md)
  : Is the IPOPT backend available?
- [`opt_classify()`](https://adamrsophie.github.io/optcompR/reference/opt_classify.md)
  : Classify an optimization model
- [`opt_standard_form()`](https://adamrsophie.github.io/optcompR/reference/opt_standard_form.md)
  : Standard form of a model
- [`opt_write_model()`](https://adamrsophie.github.io/optcompR/reference/opt_write_model.md)
  : Write a model to an MPS or LP file

## Results

- [`opt_status()`](https://adamrsophie.github.io/optcompR/reference/opt_status.md)
  [`opt_primal_status()`](https://adamrsophie.github.io/optcompR/reference/opt_status.md)
  [`opt_dual_status()`](https://adamrsophie.github.io/optcompR/reference/opt_status.md)
  [`opt_has_solution()`](https://adamrsophie.github.io/optcompR/reference/opt_status.md)
  : Solver status of a result
- [`opt_value()`](https://adamrsophie.github.io/optcompR/reference/opt_value.md)
  : Values of expressions at a solution
- [`opt_objective_value()`](https://adamrsophie.github.io/optcompR/reference/opt_objective_value.md)
  [`opt_objective_bound()`](https://adamrsophie.github.io/optcompR/reference/opt_objective_value.md)
  [`opt_gap()`](https://adamrsophie.github.io/optcompR/reference/opt_objective_value.md)
  : Objective value of a result
- [`opt_dual()`](https://adamrsophie.github.io/optcompR/reference/opt_dual.md)
  [`opt_reduced_cost()`](https://adamrsophie.github.io/optcompR/reference/opt_dual.md)
  : Dual values and reduced costs
- [`opt_slack()`](https://adamrsophie.github.io/optcompR/reference/opt_slack.md)
  : Constraint values, slacks and activity
- [`opt_solution()`](https://adamrsophie.github.io/optcompR/reference/opt_solution.md)
  : Solution table

## Diagnostics, derivatives and sensitivity

- [`opt_check()`](https://adamrsophie.github.io/optcompR/reference/opt_check.md)
  : Independent feasibility and optimality diagnostics
- [`opt_gradient()`](https://adamrsophie.github.io/optcompR/reference/opt_gradient.md)
  [`opt_jacobian()`](https://adamrsophie.github.io/optcompR/reference/opt_gradient.md)
  [`opt_hessian()`](https://adamrsophie.github.io/optcompR/reference/opt_gradient.md)
  : Exact derivatives of the objective and constraints
- [`opt_solution_jacobian()`](https://adamrsophie.github.io/optcompR/reference/opt_solution_jacobian.md)
  : Derivatives of the optimal solution with respect to parameters
- [`opt_sensitivity()`](https://adamrsophie.github.io/optcompR/reference/opt_sensitivity.md)
  : Linear programming sensitivity report
- [`opt_fixed_relaxation()`](https://adamrsophie.github.io/optcompR/reference/opt_fixed_relaxation.md)
  : Continuous problem with integer variables fixed
- [`opt_scenarios()`](https://adamrsophie.github.io/optcompR/reference/opt_scenarios.md)
  : Solve a model over parameter scenarios

## Multiobjective optimization

- [`opt_weighted()`](https://adamrsophie.github.io/optcompR/reference/opt_weighted.md)
  : Weighted-sum scalarization
- [`opt_lexicographic()`](https://adamrsophie.github.io/optcompR/reference/opt_lexicographic.md)
  : Lexicographic optimization
- [`opt_epsilon_constraint()`](https://adamrsophie.github.io/optcompR/reference/opt_epsilon_constraint.md)
  : Epsilon-constraint method for two objectives

## Integrations

- [`opt_z3_check()`](https://adamrsophie.github.io/optcompR/reference/opt_z3_check.md)
  : Exact feasibility check with Z3
- [`opt_z3_certify()`](https://adamrsophie.github.io/optcompR/reference/opt_z3_certify.md)
  : Certify a MILP or LP optimum with Z3
- [`opt_from_sym()`](https://adamrsophie.github.io/optcompR/reference/opt_from_sym.md)
  : Convert a symcompR expression to an optimization expression
- [`opt_to_sym()`](https://adamrsophie.github.io/optcompR/reference/opt_to_sym.md)
  : Convert an optimization expression to a symcompR expression
- [`opt_ode_loss()`](https://adamrsophie.github.io/optcompR/reference/opt_ode_loss.md)
  : Least squares loss of an ODE model from modelcompR

## Package

- [`optcompR`](https://adamrsophie.github.io/optcompR/reference/optcompR-package.md)
  [`optcompR-package`](https://adamrsophie.github.io/optcompR/reference/optcompR-package.md)
  : optcompR: Unified Mathematical Optimization and Computational
  Modeling in R
