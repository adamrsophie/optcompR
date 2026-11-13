#' Solve an optimization model
#'
#' Lowers the model to its standard form, classifies it, selects a solver
#' that supports the class, and returns a standardized result. With
#' `solver = "auto"` the first installed backend (in registry priority order)
#' whose R binding supports the model class is used; for nonconvex classes
#' backends with global guarantees come first. A model is never transformed
#' into an approximation to fit a solver: if no installed backend supports
#' it, an error lists what is missing.
#'
#' @param model An [opt_model()].
#' @param solver Backend name from [opt_solvers()], or `"auto"`.
#' @param time_limit Wall-clock limit in seconds, for solvers that accept it.
#' @param mip_gap Relative optimality gap for mixed-integer solvers. `NULL`
#'   keeps the solver default.
#' @param verbose Print the solver log.
#' @param warm_start Pass start values set with [opt_set_start()] or
#'   `opt_variable(start = )` to solvers that accept them.
#' @param options Named list of backend options passed through unchanged
#'   (for HiGHS, any option of [highs::highs_control()]; HiGHS defaults to
#'   one thread, and a solve with a different `threads` value than the
#'   previous one resets the HiGHS global scheduler first).
#' @return An `opt_result`. See [opt_status()], [opt_value()],
#'   [opt_objective_value()], [opt_dual()] and [opt_check()].
#' @examples
#' m <- opt_model("max")
#' x <- opt_variable(m, "x", lower = 0, upper = 2)
#' y <- opt_variable(m, "y", lower = 0, upper = 3)
#' opt_constraint(m, x + y <= 4)
#' opt_objective(m, 3 * x + 2 * y)
#' r <- opt_solve(m)
#' r
#' opt_value(r, x)
#' @export
opt_solve <- function(model, solver = "auto", time_limit = Inf, mip_gap = NULL,
                      verbose = FALSE, warm_start = TRUE, options = list()) {
  check_model(model)
  if (!is.numeric(time_limit) || length(time_limit) != 1L || time_limit <= 0) {
    stop("`time_limit` must be a positive number of seconds.", call. = FALSE)
  }
  if (!is.list(options)) stop("`options` must be a named list.", call. = FALSE)
  if (n_symbols(model) == 0L || !any(model$sym_kind == "var")) {
    stop("The model has no decision variables.", call. = FALSE)
  }
  t0 <- proc.time()[["elapsed"]]
  sf <- opt_standard_form(model)
  cls <- opt_classify(sf)
  t1 <- proc.time()[["elapsed"]]
  entry <- select_solver(cls, solver)
  control <- list(time_limit = time_limit, mip_gap = mip_gap,
                  verbose = verbose, warm_start = warm_start,
                  options = options)
  warnings <- character()
  res <- withCallingHandlers(
    entry$solve(sf, cls, control),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  t2 <- proc.time()[["elapsed"]]
  build_result(model, sf, cls, entry, res, control, warnings,
               timing = c(build = t1 - t0, solve = t2 - t1))
}

build_result <- function(model, sf, cls, entry, res, control, warnings,
                         timing) {
  z <- sf$parameter_values
  z[!is.finite(z)] <- 0
  has_x <- !anyNA(res$x)
  z[sf$var_sym] <- if (has_x) res$x else NA_real_
  names(res$x) <- sf$col_names
  obj <- NA_real_
  if (has_x) {
    obj <- if (is.null(model$objective)) 0 else
      eval_expr(model$objective, z)
  }
  if (!is.null(res$row_dual)) names(res$row_dual) <- sf$row_names
  if (!is.null(res$col_dual)) names(res$col_dual) <- sf$col_names
  term <- res$termination
  if (term == "optimal" && isTRUE(res$local)) term <- "locally_optimal"
  fd <- unique(unlist(lapply(c(list(model$objective),
                               lapply(model$cons, `[[`, "expr")),
                             function(e) if (is.null(e)) NULL else expr_fd_nodes(e))))
  if (length(fd)) {
    warnings <- c(warnings, sprintf(paste(
      "Gradients of %s are finite-difference approximations; optimality",
      "and multipliers hold only to their accuracy."),
      paste(sQuote(fd, FALSE), collapse = ", ")))
  }
  structure(list(
    model = model, solver = entry$name, solver_version = entry$version(),
    termination = term, primal_status = res$primal_status,
    dual_status = res$dual_status, objective_value = obj,
    solver_objective = res$objective, x = res$x, z = z,
    row_dual = res$row_dual, col_dual = res$col_dual,
    gap = res$gap %||% NA_real_, bound = res$bound %||% NA_real_,
    iterations = res$iterations, message = res$message,
    classification = cls, standard_form = sf, cons = model$cons,
    objective = model$objective, sense = sf$sense,
    warm_start = isTRUE(res$warm_start), extra = res$extra,
    warnings = c(warnings, res$warnings), raw = res$raw, control = control,
    timing = timing,
    metadata = list(time = Sys.time(), R = R.version.string,
                    optcompR = as.character(utils::packageVersion("optcompR")),
                    platform = R.version$platform,
                    model_id = model$id, model_key = sf$key,
                    finite_difference = fd)
  ), class = "opt_result")
}
