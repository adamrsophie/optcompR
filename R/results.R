check_result <- function(result) {
  if (!inherits(result, "opt_result")) {
    stop("`result` must be an `opt_result` returned by opt_solve().",
         call. = FALSE)
  }
  invisible(result)
}

#' Solver status of a result
#'
#' `opt_status()` returns the termination status, one of:
#'
#' * `"optimal"`: the solver proved optimality within its tolerances
#'   (globally, for the classes the backend solves globally);
#' * `"locally_optimal"`: a local solver converged; nonconvex problems may
#'   have better solutions elsewhere;
#' * `"infeasible"`, `"unbounded"`, `"infeasible_or_unbounded"`: as
#'   certified by the solver; HiGHS may report infeasible or unbounded
#'   without deciding which;
#' * `"time_limit"`, `"iteration_limit"`, `"solution_limit"`,
#'   `"objective_limit"`, `"memory_limit"`, `"interrupted"`: the solver
#'   stopped early; a feasible point may still exist, see
#'   `opt_primal_status()`;
#' * `"numerical_error"`, `"invalid_model"`, `"unknown"`.
#'
#' A mixed-integer incumbent found before a limit is never reported as
#' optimal. `opt_primal_status()` and `opt_dual_status()` describe the
#' returned points: `"feasible_point"`, `"infeasible_point"`,
#' `"no_solution"`, `"not_available"` or `"unknown"`.
#'
#' @param result An `opt_result`.
#' @return A string.
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x", lower = 0)
#' opt_constraint(m, x <= -1)
#' opt_objective(m, x)
#' opt_status(opt_solve(m))
#' @export
opt_status <- function(result) {
  check_result(result)
  result$termination
}

#' @rdname opt_status
#' @export
opt_primal_status <- function(result) {
  check_result(result)
  result$primal_status
}

#' @rdname opt_status
#' @export
opt_dual_status <- function(result) {
  check_result(result)
  result$dual_status
}

#' @rdname opt_status
#' @export
opt_has_solution <- function(result) {
  check_result(result)
  !anyNA(result$x) && result$primal_status %in% c("feasible_point", "unknown")
}

#' Values of expressions at a solution
#'
#' Evaluates any expression of the model (variables, slices, or derived
#' expressions) at the returned primal point, using the parameter values that
#' were in force when the model was solved. Later parameter changes do not
#' affect stored results.
#'
#' @param result An `opt_result`.
#' @param x An expression of the solved model, or a variable block name.
#' @return A numeric vector, or an array for matrix expressions.
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x", dim = 2, lower = 1)
#' opt_objective(m, sum(x))
#' r <- opt_solve(m)
#' opt_value(r, x)
#' opt_value(r, 2 * x[1] + 1)
#' @export
opt_value <- function(result, x) {
  check_result(result)
  if (is.character(x)) x <- opt_get(result$model, x)
  x <- as_expr(x)
  if (!is.null(x$model) && !identical(x$model, result$model)) {
    stop("The expression belongs to a different model.", call. = FALSE)
  }
  syms <- expr_symbols(x)
  if (length(syms) && max(syms) > length(result$z)) {
    stop("The expression uses variables added after the model was solved.",
         call. = FALSE)
  }
  if (!opt_has_solution(result)) {
    warning(sprintf("The result has no primal solution (status '%s').",
                    result$termination), call. = FALSE)
  }
  v <- eval_expr(x, result$z)
  if (!is.null(x$dim)) {
    v <- array(v, dim = x$dim, dimnames = x$dimnames)
  } else if (!is.null(x$dimnames) && length(x$dimnames) == 1L) {
    names(v) <- x$dimnames[[1L]]
  }
  v
}

#' Objective value of a result
#'
#' The objective value is recomputed from the model's objective expression
#' at the returned point, including constant terms, so it does not depend on
#' how a backend handles offsets. The value reported by the solver is kept
#' as `solver = TRUE`; [opt_check()] compares the two. For mixed-integer
#' problems `opt_objective_bound()` returns the best proven bound and
#' `opt_gap()` the relative gap.
#'
#' @param result An `opt_result`.
#' @param solver Return the solver-reported value instead.
#' @return A number (`NA` when there is no primal point).
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x", lower = 2)
#' opt_objective(m, 3 * x + 10)
#' opt_objective_value(opt_solve(m))
#' @export
opt_objective_value <- function(result, solver = FALSE) {
  check_result(result)
  if (solver) result$solver_objective else result$objective_value
}

#' @rdname opt_objective_value
#' @export
opt_objective_bound <- function(result) {
  check_result(result)
  result$bound
}

#' @rdname opt_objective_value
#' @export
opt_gap <- function(result) {
  check_result(result)
  result$gap
}

block_rows <- function(result, constraint) {
  name <- if (inherits(constraint, "opt_constraint")) constraint$name else
    constraint
  check_scalar_string(name, "constraint")
  sf <- result$standard_form
  bnames <- vapply(result$cons, `[[`, "", "name")
  b <- match(name, bnames)
  if (is.na(b)) stop(sprintf("No constraint block named '%s' in this result.",
                             name), call. = FALSE)
  list(rows = which(sf$row_block == b), block = result$cons[[b]])
}

#' Dual values and reduced costs
#'
#' Duals follow one convention for every backend and both senses: the dual
#' of a constraint row is the shadow price \eqn{\partial f^* / \partial b},
#' the rate of change of the optimal objective when the active bound of the
#' row is increased. The reduced cost of a variable is the same derivative
#' for its active bound. Under this convention the Lagrangian is
#' \eqn{L(x, y, z) = f(x) - \sum_i y_i (g_i(x) - b_i) - \sum_j z_j (x_j - l_j)}
#' and stationarity reads \eqn{\nabla f(x) - J(x)^\top y - z = 0}. For a
#' minimization, rows bounded above have \eqn{y \le 0} and rows bounded below
#' \eqn{y \ge 0}; signs flip for maximization.
#'
#' Duals are defined for continuous problems. For mixed-integer problems they
#' are not available, because the optimal value is not differentiable in the
#' data in general; [opt_fixed_relaxation()] builds the continuous problem
#' with integers fixed, whose duals are a labeled local quantity.
#'
#' @param result An `opt_result`.
#' @param constraint An `opt_constraint` reference or block name.
#' @return A numeric vector named by row.
#' @examples
#' m <- opt_model("max")
#' x <- opt_variable(m, "x", lower = 0)
#' y <- opt_variable(m, "y", lower = 0)
#' cap <- opt_constraint(m, x + y <= 4)
#' opt_constraint(m, x <= 2)
#' opt_objective(m, 3 * x + 2 * y)
#' r <- opt_solve(m)
#' opt_dual(r, cap)
#' opt_reduced_cost(r, x)
#' @export
opt_dual <- function(result, constraint) {
  check_result(result)
  if (is.null(result$row_dual)) {
    stop(dual_unavailable_reason(result), call. = FALSE)
  }
  br <- block_rows(result, constraint)
  result$row_dual[br$rows]
}

#' @rdname opt_dual
#' @param x A variable expression or block name.
#' @export
opt_reduced_cost <- function(result, x) {
  check_result(result)
  if (is.null(result$col_dual)) {
    stop(dual_unavailable_reason(result), call. = FALSE)
  }
  sym <- resolve_symbols(result$model, x, "var")
  result$col_dual[result$standard_form$col_of_sym[sym]]
}

dual_unavailable_reason <- function(result) {
  if (result$classification$integer) {
    return(paste("Duals are not defined for mixed-integer problems;",
                 "see opt_fixed_relaxation() for duals of the problem with",
                 "integers fixed."))
  }
  sprintf("The solver '%s' returned no dual values (status '%s').",
          result$solver, result$termination)
}

#' Constraint values, slacks and activity
#'
#' Evaluates each row of a constraint block from the original expression at
#' the solution, independently of the standard form.
#'
#' @param result An `opt_result`.
#' @param constraint An `opt_constraint` reference or block name.
#' @param tol Absolute tolerance, scaled by the bound magnitude, for marking
#'   rows active.
#' @return A data frame with the row value, lower and upper bound, slack
#'   (distance to the nearest bound, negative when violated) and an `active`
#'   flag.
#' @examples
#' m <- opt_model("max")
#' x <- opt_variable(m, "x", lower = 0)
#' cn <- opt_constraint(m, x <= 3)
#' opt_objective(m, x)
#' opt_slack(opt_solve(m), cn)
#' @export
opt_slack <- function(result, constraint, tol = 1e-7) {
  check_result(result)
  br <- block_rows(result, constraint)
  b <- br$block
  v <- eval_expr(b$expr, result$z)
  lo <- b$lower
  up <- b$upper
  s <- pmin(v - lo, up - v)
  scale <- pmax(1, abs(ifelse(is.finite(lo), lo, 0)), abs(ifelse(is.finite(up), up, 0)))
  data.frame(row = result$standard_form$row_names[br$rows],
             value = v, lower = lo, upper = up, slack = s,
             active = abs(v - lo) <= tol * scale | abs(up - v) <= tol * scale,
             stringsAsFactors = FALSE)
}

#' Solution table
#'
#' @param result An `opt_result`.
#' @return A data frame with one row per variable: name, value, bounds,
#'   type and (for continuous problems) reduced cost.
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x", dim = 2, lower = c(1, 2))
#' opt_objective(m, sum(x))
#' opt_solution(opt_solve(m))
#' @export
opt_solution <- function(result) {
  check_result(result)
  sf <- result$standard_form
  out <- data.frame(name = sf$col_names, value = unname(result$x),
                    lower = sf$lower, upper = sf$upper, type = sf$type,
                    stringsAsFactors = FALSE)
  if (!is.null(result$col_dual)) out$reduced_cost <- unname(result$col_dual)
  out
}

#' @export
print.opt_result <- function(x, ...) {
  cat("<opt_result>\n")
  cat(sprintf("  solver:      %s %s\n", x$solver, x$solver_version))
  cat(sprintf("  problem:     %s (%s%s)\n", x$classification$class,
              x$classification$convexity,
              if (x$classification$integer) " relaxation" else ""))
  cat(sprintf("  status:      %s (primal: %s, dual: %s)\n", x$termination,
              x$primal_status, x$dual_status))
  cat(sprintf("  objective:   %s\n", fmt_num(x$objective_value)))
  if (x$classification$integer && is.finite(x$bound)) {
    cat(sprintf("  bound:       %s (gap %s)\n", fmt_num(x$bound),
                fmt_num(x$gap)))
  }
  cat(sprintf("  time:        %.3fs (build %.3fs)\n", x$timing[["solve"]],
              x$timing[["build"]]))
  if (length(x$warnings)) {
    cat("  warnings:\n")
    cat(paste0("    ", x$warnings), sep = "\n")
  }
  invisible(x)
}

#' @export
summary.opt_result <- function(object, ...) {
  out <- list(result = object, diagnostics = if (opt_has_solution(object))
    opt_check(object) else NULL, solution = opt_solution(object))
  class(out) <- "summary.opt_result"
  out
}

#' @export
print.summary.opt_result <- function(x, ...) {
  print(x$result)
  cat("\nClassification evidence:\n")
  cat(paste0("  - ", x$result$classification$evidence), sep = "\n")
  if (!is.null(x$diagnostics)) {
    cat("\n")
    print(x$diagnostics)
  }
  sol <- x$solution
  cat("\nVariables:\n")
  print(utils::head(sol, 20L), row.names = FALSE)
  if (nrow(sol) > 20L) cat(sprintf("... %d more\n", nrow(sol) - 20L))
  invisible(x)
}
