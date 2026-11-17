#' Solve an optimization model
#'
#' Lowers the model to its standard form, classifies it, selects a solver
#' that supports the class, and returns a standardized result. With
#' `solver = "auto"` the first installed backend (in registry priority order)
#' whose R binding supports the model class is used; for nonconvex classes
#' backends with global guarantees come first, and convex QPs go to the
#' Clarabel interior point solver before HiGHS when it is installed. A model is never transformed
#' into an approximation to fit a solver: if no installed backend supports
#' it, an error lists what is missing.
#'
#' @param model An [opt_model()].
#' @param solver Backend name from [opt_solvers()], or `"auto"`.
#' @param time_limit Wall-clock limit in seconds, for solvers that accept it.
#' @param mip_gap Relative optimality gap for mixed-integer solvers. `NULL`
#'   keeps the solver default. A solve that stops with a gap above zero
#'   under a gap the user asked for (here or through a native option such
#'   as HiGHS `mip_rel_gap` or SCIP `limits/gap`), or with a gap above 1e-4
#'   in any case, has status `"gap_limit"`, not `"optimal"`.
#' @param verbose Print the solver log.
#' @param warm_start Pass start values set with [opt_set_start()] or
#'   `opt_variable(start = )` to solvers that accept them, and reuse the
#'   native HiGHS solver and OSQP workspace of the previous solve of the
#'   same model: when only costs, bounds or matrix values changed (for
#'   example after [opt_set_parameter()]), they are updated in place and the
#'   solver restarts from its last basis or solution. `FALSE` builds a new
#'   native solver.
#' @param options Named list of backend options passed through unchanged
#'   (for HiGHS, any option of [highs::highs_control()]). HiGHS uses one
#'   thread by default; its thread scheduler is global to the process, so
#'   the `threads` value of the first HiGHS solve holds for the whole R
#'   session and a different value is an error.
#' @section Checking the returned point:
#' The point a backend returns is checked against the model (rows from the
#' original expressions, bounds, integrality, cones, indicators and SOS
#' constraints), each side scaled by `max(1, |bound|)`, with the backend's
#' feasibility tolerance (1e-5 by default, 1e-4 for IPOPT and 1e-2 for the
#' first-order solvers OSQP and SCS); cones, indicators and SOS
#' constraints are checked at its square root. When `options` request a
#' stricter accuracy from OSQP or SCS (`eps_abs`, `eps_rel`), the tolerance
#' becomes 100 times the strictest request, so the check matches the
#' accuracy asked for. The resolved tolerances are stored in
#' `result$validation`. A point that fails is labeled
#' `"infeasible_point"`, a claimed optimum then becomes
#' `"locally_infeasible"` (local solvers) or `"numerical_error"`, and a
#' warning gives the violation. An objective that is not finite at the
#' point (for example `log(0)`) also makes it `"infeasible_point"`. A
#' backend is treated as local for the classes it does not declare global
#' unless the model is continuous and proven convex.
#'
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
  entry <- select_solver(cls, solver, sf)
  # per model and backend scratch space, for example a basis to restart from
  if (is.null(model$cache$backend)) model$cache$backend <- new.env(parent = emptyenv())
  if (is.null(model$cache$backend[[entry$name]])) {
    model$cache$backend[[entry$name]] <- new.env(parent = emptyenv())
  }
  control <- list(time_limit = time_limit, mip_gap = mip_gap,
                  verbose = verbose, warm_start = warm_start,
                  options = options, cache = model$cache$backend[[entry$name]])
  warnings <- unsupported_controls(entry, control, sf)
  sf_solver <- entry$prepared_sf %||% reformulate_special(sf, entry)
  entry$prepared_sf <- NULL
  if (!is.null(sf_solver$reformulated)) {
    warnings <- c(warnings, sprintf(paste(
      "Indicator/SOS constraints were reformulated exactly for %s with %d",
      "binary column(s) and %d row(s)."), entry$name,
      sf_solver$reformulated$added_columns, sf_solver$reformulated$added_rows))
  }
  res <- withCallingHandlers(
    trim_result(check_backend_result(entry$solve(sf_solver, cls, control),
                                     entry, sf_solver), sf),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  t2 <- proc.time()[["elapsed"]]
  build_result(model, sf, cls, entry, res, control, warnings,
               timing = c(build = t1 - t0, solve = t2 - t1))
}

# Controls the chosen backend does not declare are ignored by it; say so.
unsupported_controls <- function(entry, control, sf) {
  f <- entry$features
  c(if (!is.null(control$mip_gap) && !isTRUE(f$mip_gap))
      sprintf("%s does not take a MIP gap; `mip_gap` was ignored.", entry$name),
    if (is.finite(control$time_limit) && !isTRUE(f$time_limit))
      sprintf("%s does not take a time limit; `time_limit` was ignored.",
              entry$name),
    if (isTRUE(control$warm_start) && any(!is.na(sf$start)) &&
        !isTRUE(f$warm_start))
      sprintf("%s does not use start values; they were ignored.", entry$name))
}

# The tolerances a returned point is validated against: the backend's
# default, tightened to the accuracy the user asked the solver for (for
# SCS and OSQP, eps_abs and eps_rel), so a strict request is not met by a
# percent-level check. Rows and bounds use `tol`, cone and other special
# constraints `special_tol`.
validation_policy <- function(entry, control) {
  tol <- entry$feasibility_tol %||% 1e-5
  source <- "solver default"
  opts <- control$options[intersect(names(control$options),
                                    entry$accuracy_options)]
  req <- suppressWarnings(as.numeric(unlist(opts)))
  req <- req[is.finite(req) & req > 0]
  if (length(req)) {
    t2 <- max(100 * max(req), 1e-8)
    if (t2 < tol) {
      tol <- t2
      source <- "requested solver accuracy"
    }
  }
  list(tol = tol, special_tol = sqrt(tol), source = source)
}

build_result <- function(model, sf, cls, entry, res, control, warnings,
                         timing) {
  z <- sf$parameter_values
  z[!is.finite(z)] <- 0
  has_x <- !anyNA(res$x)
  z[sf$var_sym] <- if (has_x) res$x else NA_real_
  names(res$x) <- sf$col_names
  obj <- NA_real_
  obj_warn <- character()
  if (has_x) {
    obj <- if (is.null(model$objective)) 0 else withCallingHandlers(
      eval_expr(model$objective, z), warning = function(w) {
        obj_warn <<- c(obj_warn, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
  }
  if (!is.null(res$row_dual)) names(res$row_dual) <- sf$row_names
  if (!is.null(res$col_dual)) names(res$col_dual) <- sf$col_names
  term <- res$termination
  # a backend is local for the classes it does not declare global, unless
  # the model is continuous and proven convex, where a local optimum is
  # global; for a mixed-integer model the convexity is the relaxation's,
  # which does not make a local search global
  local <- isTRUE(res$local) || (!cls$class %in% entry$global &&
    (cls$integer || cls$convexity != "convex"))
  # the returned point is checked against the model: a backend's
  # feasibility claim is not trusted on its own
  validation <- validation_policy(entry, control)
  if (has_x && res$primal_status %in% c("feasible_point", "unknown")) {
    tol <- validation$tol
    fs <- feasibility_summary(point_feasibility(sf, model$cons, z, tol), tol)
    if (!fs$ok || any(!is.finite(res$x))) {
      res$primal_status <- "infeasible_point"
      warnings <- c(warnings, sprintf(paste(
        "%s returned a point that violates the model (largest scaled",
        "violation %s); it is labeled infeasible_point."), entry$name,
        fmt_num(max(fs$max))))
      if (term %in% c("optimal", "locally_optimal", "almost_optimal")) {
        term <- if (local) "locally_infeasible" else "numerical_error"
      }
    } else if (res$primal_status == "unknown") {
      res$primal_status <- "feasible_point"
    }
  }
  # the objective must be defined at a usable point: log(0) or sqrt(-1)
  # is a domain failure, not an optimum; the point is kept for diagnosis
  if (has_x && !is.finite(obj) &&
      res$primal_status %in% c("feasible_point", "unknown")) {
    res$primal_status <- "infeasible_point"
    warnings <- c(warnings, sprintf(paste(
      "The objective is %s at the point %s returned%s; it is labeled",
      "infeasible_point."), format(obj), entry$name,
      if (length(obj_warn)) sprintf(" (%s)", paste(unique(obj_warn), collapse = "; ")) else ""))
    if (term %in% c("optimal", "locally_optimal", "almost_optimal")) {
      term <- "numerical_error"
    }
  } else if (length(obj_warn)) {
    warnings <- c(warnings, sprintf("Evaluating the objective warned: %s",
                                    paste(unique(obj_warn), collapse = "; ")))
  }
  if (term == "optimal" && local) term <- "locally_optimal"
  # a stop at a relative gap is not a proof of optimality: a mixed-integer
  # result with a gap above 1e-9 under a gap the user requested (through
  # `mip_gap` or a native gap option), or above 1e-4 (looser than the
  # solvers' defaults) in any case, has status gap_limit
  if (term == "optimal" && cls$integer && is.finite(res$gap %||% NA_real_)) {
    requested <- !is.null(control$mip_gap) || any(names(control$options) %in%
      c("mip_rel_gap", "mip_abs_gap", "limits/gap", "limits/absgap"))
    if (res$gap > 1e-4 || (requested && res$gap > 1e-9)) term <- "gap_limit"
  }
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
    model = model, model_key = paste(model$structure_version,
                                     model$parameter_version, sep = ":"),
    solver = entry$name, solver_version = entry$version(),
    termination = term, primal_status = res$primal_status,
    dual_status = res$dual_status, objective_value = obj,
    solver_objective = res$objective, x = res$x, z = z,
    row_dual = res$row_dual, col_dual = res$col_dual,
    gap = res$gap %||% NA_real_, bound = res$bound %||% NA_real_,
    iterations = res$iterations, message = res$message,
    classification = cls, standard_form = sf, cons = model$cons,
    objective = model$objective, special = model$special, sense = sf$sense,
    warm_start = isTRUE(res$warm_start), extra = res$extra,
    validation = validation,
    warnings = c(warnings, res$warnings), raw = res$raw, control = control,
    timing = timing,
    metadata = list(time = Sys.time(), R = R.version.string,
                    optcompR = as.character(utils::packageVersion("optcompR")),
                    platform = R.version$platform,
                    model_id = model$id, model_key = sf$key,
                    finite_difference = fd)
  ), class = "opt_result")
}
