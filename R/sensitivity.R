#' Linear programming sensitivity report
#'
#' For an optimal linear program solved by HiGHS, reports the ranges over
#' which objective coefficients and row bounds can move while the optimal
#' basis stays optimal, together with the objective value at the ends of
#' each range. Ranges are local to the returned basis; at degenerate optima
#' different bases give different ranges, and outside a range the optimal
#' value is piecewise linear rather than following the dual.
#'
#' Row bounds are reported in standard-form terms (`row_lower <= a'x <=
#' row_upper`, constants moved to the bounds), so for `a'x + c <= b` the
#' bound is `b - c`.
#'
#' @param result An `opt_result` for an LP solved with `solver = "highs"`.
#' @return A list of two data frames: `costs` (variable, coefficient, lowest
#'   and highest coefficient keeping the basis optimal, and the objective
#'   at those values) and `rows` (row, active bound, dual, lowest and highest
#'   bound value and the objective there).
#' @examples
#' m <- opt_model("max")
#' x <- opt_variable(m, "x", lower = 0)
#' y <- opt_variable(m, "y", lower = 0)
#' opt_constraint(m, x + y <= 4)
#' opt_constraint(m, x <= 2)
#' opt_objective(m, 3 * x + 2 * y)
#' opt_sensitivity(opt_solve(m))
#' @export
opt_sensitivity <- function(result) {
  check_result(result)
  rg <- result$extra$ranging
  if (is.null(rg) || !isTRUE(rg$valid)) {
    stop(paste("Sensitivity ranging is available for optimal linear programs",
               "solved with HiGHS."), call. = FALSE)
  }
  sf <- result$standard_form
  s <- result$extra$sign
  n <- length(sf$c)
  k <- nrow(sf$A)
  take <- function(part, len) lapply(part, function(v) v[seq_len(len)])
  cu <- take(rg$col_cost_up, n)
  cd <- take(rg$col_cost_dn, n)
  # solver works with s * c; map back to user coefficients
  if (s > 0) {
    lo_c <- cd$value; lo_o <- cd$objective
    hi_c <- cu$value; hi_o <- cu$objective
  } else {
    lo_c <- -cu$value; lo_o <- -cu$objective
    hi_c <- -cd$value; hi_o <- -cd$objective
  }
  costs <- data.frame(variable = sf$col_names, coefficient = sf$c,
                      lower = lo_c, upper = hi_c,
                      objective_at_lower = lo_o, objective_at_upper = hi_o,
                      stringsAsFactors = FALSE)
  rows <- data.frame(row = character(), stringsAsFactors = FALSE)
  if (k > 0L) {
    ru <- take(rg$row_bound_up, k)
    rd <- take(rg$row_bound_dn, k)
    rows <- data.frame(row = sf$row_names, row_lower = sf$row_lower,
                       row_upper = sf$row_upper,
                       dual = unname(result$row_dual),
                       bound_lower = rd$value, bound_upper = ru$value,
                       objective_at_lower = s * rd$objective,
                       objective_at_upper = s * ru$objective,
                       stringsAsFactors = FALSE)
  }
  list(costs = costs, rows = rows)
}

#' Continuous problem with integer variables fixed
#'
#' Fixes every integer and binary variable of a mixed-integer result at its
#' returned value and solves the remaining continuous problem. Its duals
#' describe the sensitivity of the continuous part with the discrete choices
#' held fixed. They are not duals of the mixed-integer problem, whose optimal
#' value is in general not differentiable in the data; the returned result
#' is labeled with `fixed_relaxation = TRUE`.
#'
#' @param result An `opt_result` of a mixed-integer model with a solution.
#' @param ... Further arguments for [opt_solve()].
#' @return An `opt_result` for the continuous problem. The original model is
#'   not changed.
#' @examples
#' m <- opt_model("max")
#' x <- opt_variable(m, "x", type = "integer", lower = 0, upper = 3)
#' y <- opt_variable(m, "y", lower = 0)
#' cap <- opt_constraint(m, x + y <= 4.5)
#' opt_objective(m, 3 * x + 2 * y)
#' r <- opt_solve(m)
#' fr <- opt_fixed_relaxation(r)
#' opt_dual(fr, cap)
#' @export
opt_fixed_relaxation <- function(result, ...) {
  check_result(result)
  if (!result$classification$integer) {
    stop("The result is not from a mixed-integer model.", call. = FALSE)
  }
  if (!opt_has_solution(result)) {
    stop("The result has no integer solution to fix.", call. = FALSE)
  }
  model <- result$model
  vs <- result$standard_form$var_sym
  int <- vs[model$type[vs] %in% c("integer", "binary")]
  saved <- list(lower = model$lower[int], upper = model$upper[int],
                type = model$type[int])
  vals <- round(result$z[int])
  on.exit({
    model$lower[int] <- saved$lower
    model$upper[int] <- saved$upper
    model$type[int] <- saved$type
    touch_structure(model)
  }, add = TRUE)
  model$lower[int] <- vals
  model$upper[int] <- vals
  model$type[int] <- "continuous"
  touch_structure(model)
  out <- opt_solve(model, ...)
  out$fixed_relaxation <- TRUE
  out$warnings <- c(out$warnings, paste(
    "Continuous problem with integer variables fixed at the incumbent;",
    "duals are conditional on the discrete choices."))
  out
}

#' Solve a model over parameter scenarios
#'
#' Sets a parameter to each value in turn, solves, and collects the status
#' and objective. The parameter is restored afterwards. Each solve rebuilds
#' the numerical data from the cached structure.
#'
#' @param model An [opt_model()].
#' @param parameter A parameter expression or name.
#' @param values A list of values or a numeric vector (one scenario per
#'   element for scalar parameters).
#' @param ... Further arguments for [opt_solve()].
#' @param keep Keep the full results in the `results` attribute.
#' @return A data frame with one row per scenario: index, status and
#'   objective value.
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x", lower = 0)
#' d <- opt_parameter(m, "d", value = 1)
#' opt_constraint(m, x >= d)
#' opt_objective(m, 2 * x)
#' opt_scenarios(m, d, c(1, 2, 3))
#' @export
opt_scenarios <- function(model, parameter, values, ..., keep = FALSE) {
  check_model(model)
  sym <- resolve_symbols(model, parameter, "par")
  old <- model$value[sym]
  on.exit(opt_set_parameter(model, sym_expr_for(model, sym), old), add = TRUE)
  if (!is.list(values)) values <- as.list(values)
  results <- vector("list", length(values))
  for (k in seq_along(values)) {
    model$value[sym] <- recycle_to(as.numeric(values[[k]]), length(sym),
                                   "value")
    touch_parameters(model)
    results[[k]] <- opt_solve(model, ...)
  }
  out <- data.frame(
    scenario = seq_along(values),
    value = vapply(values, function(v) paste(fmt_num(v), collapse = ","), ""),
    status = vapply(results, opt_status, ""),
    objective = vapply(results, opt_objective_value, 0),
    stringsAsFactors = FALSE)
  if (keep) attr(out, "results") <- results
  out
}

sym_expr_for <- function(model, sym) symbol_expr(model, sym, NULL)
