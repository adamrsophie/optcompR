# Multiobjective optimization by scalarization. Every function restores the
# model's objective and removes temporary constraints on exit.

with_objective <- function(model, expr, code) {
  old <- model$objective
  on.exit({
    model$objective <- old
    touch_structure(model)
  }, add = TRUE)
  opt_objective(model, expr)
  force(code)
}

check_objectives <- function(model, objectives) {
  if (!is.list(objectives) || length(objectives) < 2L) {
    stop("`objectives` must be a list of at least two scalar expressions.",
         call. = FALSE)
  }
  lapply(objectives, function(o) {
    o <- as_expr(o)
    if (o$len != 1L) stop("Each objective must be scalar.", call. = FALSE)
    if (!is.null(o$model) && !identical(o$model, model)) {
      stop("An objective belongs to a different model.", call. = FALSE)
    }
    o
  })
}

#' Weighted-sum scalarization
#'
#' Solves the model with the objective \eqn{\sum_k w_k f_k(x)} in the model's
#' sense. With positive weights every solution is Pareto optimal. Weighted
#' sums cannot reach Pareto points on nonconvex parts of the front (and in
#' mixed-integer problems they miss unsupported efficient points); use
#' [opt_epsilon_constraint()] for those.
#'
#' @param model An [opt_model()].
#' @param objectives A list of scalar expressions.
#' @param weights Nonnegative weights, one per objective.
#' @param ... Further arguments for [opt_solve()].
#' @return An `opt_result`, with the individual objective values in
#'   `objective_values`.
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x", lower = 0, upper = 4)
#' y <- opt_variable(m, "y", lower = 0, upper = 4)
#' opt_constraint(m, x + y >= 4)
#' r <- opt_weighted(m, list(cost = 2 * x + y, risk = x + 3 * y), c(1, 1))
#' r$objective_values
#' @export
opt_weighted <- function(model, objectives, weights, ...) {
  check_model(model)
  objs <- check_objectives(model, objectives)
  check_numeric(weights, "weights")
  if (length(weights) != length(objs) || any(weights < 0)) {
    stop("`weights` must be nonnegative, one per objective.", call. = FALSE)
  }
  total <- Reduce(expr_add, Map(function(o, w) expr_scale(o, w), objs, weights))
  r <- with_objective(model, total, opt_solve(model, ...))
  r$objective_values <- objective_values(r, objs, objectives)
  r
}

objective_values <- function(r, objs, original) {
  vals <- if (opt_has_solution(r)) vapply(objs, function(o) eval_expr(o, r$z), 0)
    else rep(NA_real_, length(objs))
  stats::setNames(vals, names(original))
}

#' Lexicographic optimization
#'
#' Optimizes the objectives in priority order. After objective \eqn{k} is
#' optimized with value \eqn{f_k^*}, the constraint
#' \eqn{f_k(x) \le f_k^* + \mathrm{tol} \max(1, |f_k^*|)} (or \eqn{\ge} for
#' maximization) is added before the next objective. The relative tolerance
#' keeps later stages feasible despite solver tolerances; it also means
#' later stages may trade a tiny amount of earlier objectives.
#'
#' @inheritParams opt_weighted
#' @param tol Relative tolerance on earlier objectives.
#' @return A list of `opt_result` objects, one per stage; the last one is
#'   the lexicographic solution. Objective values of all stages are in its
#'   `objective_values`.
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x", lower = 0, upper = 4)
#' y <- opt_variable(m, "y", lower = 0, upper = 4)
#' opt_constraint(m, x + y >= 4)
#' st <- opt_lexicographic(m, list(cost = x + y, risk = x))
#' opt_value(st[[2]], c(x, y))
#' @export
opt_lexicographic <- function(model, objectives, tol = 1e-7, ...) {
  check_model(model)
  objs <- check_objectives(model, objectives)
  added <- character()
  on.exit(for (nm in added) opt_delete_constraint(model, nm), add = TRUE)
  stages <- vector("list", length(objs))
  for (k in seq_along(objs)) {
    r <- with_objective(model, objs[[k]], opt_solve(model, ...))
    stages[[k]] <- r
    if (!opt_status(r) %in% c("optimal", "locally_optimal")) {
      warning(sprintf("Stage %d ended with status '%s'; stopping.", k,
                      opt_status(r)), call. = FALSE)
      return(stages[seq_len(k)])
    }
    if (k < length(objs)) {
      v <- opt_objective_value(r)
      slack <- tol * max(1, abs(v))
      nm <- paste0(".lexicographic_", k, "_", model$id)
      rel <- if (model$sense == "min") new_relation(objs[[k]], v + slack, "<=")
        else new_relation(objs[[k]], v - slack, ">=")
      opt_constraint(model, rel, name = nm)
      added <- c(added, nm)
    }
  }
  stages[[length(objs)]]$objective_values <-
    objective_values(stages[[length(objs)]], objs, objectives)
  stages
}

#' Epsilon-constraint method for two objectives
#'
#' Traces a Pareto front by optimizing `primary` subject to
#' `secondary <= epsilon` (for minimization; `>=` for maximization) over a
#' grid of `epsilon` values. Unlike weighted sums, this reaches points on
#' nonconvex parts of the front. Points can be weakly Pareto optimal; a
#' second lexicographic stage would remove dominated ties.
#'
#' @inheritParams opt_weighted
#' @param primary,secondary Scalar expressions.
#' @param epsilon Numeric vector of bounds on `secondary`.
#' @return A data frame with `epsilon`, status, the two objective values,
#'   with the results in the `results` attribute.
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x", lower = 0, upper = 4)
#' y <- opt_variable(m, "y", lower = 0, upper = 4)
#' opt_constraint(m, x + y >= 4)
#' opt_epsilon_constraint(m, 2 * x + y, x + 3 * y, epsilon = c(4, 8, 12))
#' @export
opt_epsilon_constraint <- function(model, primary, secondary, epsilon, ...) {
  check_model(model)
  objs <- check_objectives(model, list(primary, secondary))
  check_numeric(epsilon, "epsilon")
  nm <- paste0(".epsilon_", model$id)
  res <- vector("list", length(epsilon))
  for (k in seq_along(epsilon)) {
    rel <- if (model$sense == "min") new_relation(objs[[2L]], epsilon[k], "<=")
      else new_relation(objs[[2L]], epsilon[k], ">=")
    opt_constraint(model, rel, name = nm)
    res[[k]] <- tryCatch(with_objective(model, objs[[1L]], opt_solve(model, ...)),
                         finally = opt_delete_constraint(model, nm))
  }
  vals <- t(vapply(res, function(r) objective_values(r, objs, list(1, 2)),
                   numeric(2)))
  out <- data.frame(epsilon = epsilon,
                    status = vapply(res, opt_status, ""),
                    primary = vals[, 1L], secondary = vals[, 2L],
                    stringsAsFactors = FALSE)
  attr(out, "results") <- res
  out
}
