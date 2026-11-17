# Capability-aware solver registry.

.solver_registry <- new.env(parent = emptyenv())

# Register a backend. `classes` lists the problem classes the R binding
# supports; `convex_only` lists classes accepted only when the
# classification proves (or numerically supports) convexity.
register_solver <- function(name, package, description, classes,
                            convex_only = character(), global = character(),
                            features = list(), solve, version, priority = 50L) {
  .solver_registry[[name]] <- list(
    name = name, package = package, description = description,
    classes = classes, convex_only = convex_only, global = global,
    features = utils::modifyList(list(
      duals = FALSE, ranging = FALSE, warm_start = FALSE,
      time_limit = FALSE, mip_gap = FALSE, iteration_limit = FALSE,
      cones = character(), special = character()),
      features),
    solve = solve, version = version, priority = priority)
}

solver_available <- function(entry) {
  all(vapply(entry$package, requireNamespace, logical(1), quietly = TRUE))
}

#' Available optimization solvers
#'
#' Lists the registered backends, whether their R packages are installed,
#' and the problem classes each one accepts. Capabilities describe what the
#' R binding exposes, which can be less than the native solver supports.
#'
#' @return A data frame with one row per backend.
#' @examples
#' opt_solvers()
#' @export
opt_solvers <- function() {
  nm <- names(.solver_registry)
  ent <- lapply(nm, function(n) .solver_registry[[n]])
  ord <- order(vapply(ent, `[[`, 0L, "priority"))
  ent <- ent[ord]
  data.frame(
    solver = vapply(ent, `[[`, "", "name"),
    package = vapply(ent, function(e) paste(e$package, collapse = ", "), ""),
    available = vapply(ent, solver_available, logical(1)),
    classes = vapply(ent, function(e) paste(e$classes, collapse = " "), ""),
    duals = vapply(ent, function(e) e$features$duals, logical(1)),
    description = vapply(ent, `[[`, "", "description"),
    stringsAsFactors = FALSE
  )
}

#' Capabilities of one solver backend
#'
#' @param solver Backend name, see [opt_solvers()].
#' @return A list with the problem classes, the classes accepted only when
#'   convex, the classes solved to global optimality, feature flags (duals,
#'   ranging, warm starts, time limits, MIP gaps, iteration limits),
#'   availability and the installed version.
#' @examples
#' opt_solver_capabilities("highs")
#' @export
opt_solver_capabilities <- function(solver) {
  e <- solver_entry(solver)
  avail <- solver_available(e)
  list(name = e$name, package = e$package, available = avail,
       version = if (avail) e$version() else NA_character_,
       classes = e$classes, convex_only = e$convex_only,
       global = e$global, features = e$features,
       description = e$description)
}

solver_entry <- function(solver) {
  check_scalar_string(solver, "solver")
  e <- .solver_registry[[solver]]
  if (is.null(e)) {
    stop(sprintf("Unknown solver '%s'. Registered solvers: %s.", solver,
                 paste(sort(names(.solver_registry)), collapse = ", ")),
         call. = FALSE)
  }
  e
}

# Can backend `e` solve a model of classification `cls`? Returns NULL when
# it can, or a reason string.
solver_refusal <- function(e, cls) {
  if (!cls$class %in% e$classes) {
    return(sprintf("%s does not support %s problems", e$name, cls$class))
  }
  bad_cone <- setdiff(cls$cone_types %||% character(),
                      e$features$cones %||% character())
  if (length(bad_cone)) {
    return(sprintf("%s does not support %s cone constraints", e$name,
                   paste(bad_cone, collapse = ", ")))
  }
  if (cls$class %in% e$convex_only && cls$convexity != "convex") {
    return(sprintf("%s requires a convex %s, and this model is %s",
                   e$name, cls$class,
                   if (cls$convexity == "nonconvex") "nonconvex" else
                     "not proven convex"))
  }
  NULL
}

select_solver <- function(cls, solver, sf = NULL) {
  if (!identical(solver, "auto")) {
    e <- solver_entry(solver)
    if (!solver_available(e)) {
      stop(sprintf("Solver '%s' needs the package(s) %s, which are not installed.",
                   solver, paste(e$package, collapse = ", ")), call. = FALSE)
    }
    why <- solver_refusal(e, cls)
    if (!is.null(why)) stop(why, ".", call. = FALSE)
    return(e)
  }
  ent <- as.list(.solver_registry)
  ent <- ent[order(vapply(ent, `[[`, 0L, "priority"))]
  reasons <- character()
  # prefer solvers that give global guarantees for nonconvex classes
  if (cls$convexity != "convex") {
    glob <- Filter(function(e) cls$class %in% e$global, ent)
    ent <- c(glob, Filter(function(e) !cls$class %in% e$global, ent))
  }
  for (e in ent) {
    why <- solver_refusal(e, cls)
    if (!is.null(why)) {
      reasons <- c(reasons, why)
      next
    }
    if (!solver_available(e)) {
      reasons <- c(reasons, sprintf("%s supports %s but package %s is not installed",
                                    e$name, cls$class,
                                    paste(e$package, collapse = ", ")))
      next
    }
    # indicator and SOS constraints must be native or reformulable
    if (!is.null(sf) && (length(sf$indicators) || length(sf$sos))) {
      why <- tryCatch({
        reformulate_special(sf, e)
        NULL
      }, error = function(err) conditionMessage(err))
      if (!is.null(why)) {
        reasons <- c(reasons, why)
        next
      }
    }
    return(e)
  }
  stop(sprintf("No installed solver supports this %s (%s). %s.", cls$class,
               cls$convexity, paste(unique(reasons), collapse = "; ")),
       call. = FALSE)
}

#' Register a solver backend
#'
#' Adds a backend to the registry used by [opt_solve()] and
#' [opt_solvers()], so packages and users can connect other solvers without
#' changing optcompR. Registering an existing name replaces that backend.
#'
#' @section The solve function:
#' `solve(sf, cls, control)` receives the [opt_standard_form()] `sf` of the
#' model (parameters substituted; for indicator and SOS constraints that the
#' backend does not list in `special`, already reformulated with binary
#' columns and rows), the [opt_classify()] result `cls`, and `control`, a
#' list with `time_limit` (seconds, `Inf` for none), `mip_gap` (`NULL` for
#' the solver default), `verbose`, `warm_start`, `options` (the
#' `options` list given to [opt_solve()]) and `cache`, an environment kept
#' per model and backend in which a backend may store state, such as a
#' basis, for the next solve of the same model. `sf$structure_key` changes
#' whenever the model structure does. The objective is
#' `0.5 x'Qx + c'x + offset` in the sense `sf$sense`, or, when
#' `sf$objective_kind` is `"nonlinear"`, the lowered expression
#' `sf$objective_lowered`, which [opt_gradient()] and the other derivative
#' functions evaluate. It returns a list with
#'
#' * `termination`: one of the statuses documented in [opt_status()];
#' * `primal_status` and `dual_status`: `"feasible_point"`,
#'   `"infeasible_point"`, `"no_solution"`, `"not_available"` or
#'   `"unknown"`;
#' * `x`: the column values (length `length(sf$c)`), `NA` when there is no
#'   point;
#' * `objective`: the objective value including `sf$offset`, or `NA`;
#' * optionally `row_dual` and `col_dual` (shadow prices
#'   \eqn{\partial f^* / \partial b} of the rows and of the column bounds,
#'   in the user's sense), `gap`, `bound`, `iterations` (a named numeric
#'   vector), `message`, `local` (`TRUE` for local solvers, which turns
#'   `"optimal"` into `"locally_optimal"`), `warm_start`, `warnings` (a
#'   character vector), `extra` (a list) and `raw` (the solver's own
#'   output).
#'
#' The result is checked: a missing or unknown status, or an `x` of the
#' wrong length, is an error naming the backend.
#'
#' @param name Backend name used in `opt_solve(solver = )`.
#' @param solve The solve function described below.
#' @param classes Problem classes the backend accepts (see
#'   [opt_classify()]).
#' @param package R packages the backend needs; the backend is available
#'   when all of them are installed. `character()` for none.
#' @param description One line shown by [opt_solvers()].
#' @param convex_only Classes accepted only when convexity is proven.
#' @param global Classes the backend solves to global optimality; for other
#'   classes it is treated as a local solver.
#' @param features Named list of capabilities: `duals`, `ranging`,
#'   `warm_start`, `time_limit`, `mip_gap`, `iteration_limit` (logical),
#'   `cones` (cone types such as `"soc"`, `"exp"`, `"pow"`, `"psd"`) and
#'   `special` (`"indicator"`, `"sos"`: handled natively, otherwise
#'   reformulated).
#' @param version A function returning the solver version string, or a
#'   string.
#' @param priority Order in automatic selection; lower is tried first. The
#'   built-in backends use 10 to 90.
#' @return The backend name, invisibly. `opt_unregister_solver()` returns
#'   `TRUE` when a backend was removed.
#' @examples
#' # a backend that hands LPs to HiGHS through the highs package
#' my_lp <- function(sf, cls, control) {
#'   s <- if (sf$sense == "max") -1 else 1
#'   r <- highs::highs_solve(L = s * sf$c, lower = sf$lower,
#'                           upper = sf$upper, A = sf$A,
#'                           lhs = sf$row_lower, rhs = sf$row_upper)
#'   ok <- r$status == 7L
#'   list(termination = if (ok) "optimal" else "unknown",
#'        primal_status = if (ok) "feasible_point" else "no_solution",
#'        dual_status = "not_available",
#'        x = if (ok) r$primal_solution else rep(NA_real_, length(sf$c)),
#'        objective = if (ok) s * r$objective_value + sf$offset else NA_real_)
#' }
#' opt_register_solver("my_lp", my_lp, classes = "LP", package = "highs",
#'                     description = "HiGHS LP through a user backend")
#' m <- opt_model("max")
#' x <- opt_variable(m, "x", lower = 0, upper = 3)
#' opt_objective(m, 2 * x)
#' opt_objective_value(opt_solve(m, solver = "my_lp"))
#' opt_unregister_solver("my_lp")
#' @export
opt_register_solver <- function(name, solve, classes, package = character(),
                                description = name, convex_only = character(),
                                global = classes, features = list(),
                                version = NA_character_, priority = 50L) {
  check_scalar_string(name, "name")
  if (!is.function(solve)) stop("`solve` must be a function.", call. = FALSE)
  fa <- names(formals(solve))
  if (length(fa) < 3L && !"..." %in% fa) {
    stop("`solve` must take the arguments (sf, cls, control).", call. = FALSE)
  }
  known <- c("LP", "MILP", "QP", "MIQP", "QCQP", "MIQCQP", "SOCP", "MISOCP",
             "CONIC", "MICONIC", "NLP", "MINLP")
  if (!is.character(classes) || !length(classes) || !all(classes %in% known)) {
    stop(sprintf("`classes` must be a subset of %s.",
                 paste(known, collapse = ", ")), call. = FALSE)
  }
  if (!is.character(package)) stop("`package` must be a character vector.",
                                   call. = FALSE)
  check_scalar_string(description, "description")
  if (!is.list(features) || (length(features) && is.null(names(features)))) {
    stop("`features` must be a named list.", call. = FALSE)
  }
  bad <- setdiff(names(features), c("duals", "ranging", "warm_start",
                                    "time_limit", "mip_gap",
                                    "iteration_limit", "cones", "special"))
  if (length(bad)) stop(sprintf("Unknown feature(s): %s.",
                                paste(bad, collapse = ", ")), call. = FALSE)
  ver <- if (is.function(version)) version else {
    v <- as.character(version)
    function() v
  }
  if (!is.numeric(priority) || length(priority) != 1L) {
    stop("`priority` must be a number.", call. = FALSE)
  }
  register_solver(name, package, description, classes,
                  convex_only = intersect(convex_only, classes),
                  global = intersect(global, classes), features = features,
                  solve = solve, version = ver,
                  priority = as.integer(priority))
  invisible(name)
}

#' @rdname opt_register_solver
#' @export
opt_unregister_solver <- function(name) {
  check_scalar_string(name, "name")
  if (is.null(.solver_registry[[name]])) return(invisible(FALSE))
  rm(list = name, envir = .solver_registry)
  invisible(TRUE)
}

termination_values <- c(
  "optimal", "locally_optimal", "almost_optimal", "locally_infeasible",
  "infeasible", "unbounded", "infeasible_or_unbounded", "time_limit",
  "iteration_limit", "node_limit", "solution_limit", "objective_limit",
  "memory_limit", "interrupted", "gap_limit", "numerical_error",
  "invalid_model", "unknown")

point_values <- c("feasible_point", "infeasible_point", "no_solution",
                  "not_available", "unknown")

# Check what a backend returned before building the result.
check_backend_result <- function(res, entry, sf) {
  fail <- function(msg) {
    stop(sprintf("Backend '%s' returned %s.", entry$name, msg), call. = FALSE)
  }
  if (!is.list(res)) fail("something other than a list")
  if (!is.character(res$termination) || length(res$termination) != 1L ||
      !res$termination %in% termination_values) {
    fail(sprintf("an unknown termination status (%s)",
                 paste(format(res$termination), collapse = " ")))
  }
  for (f in c("primal_status", "dual_status")) {
    v <- res[[f]] %||% "unknown"
    if (!is.character(v) || length(v) != 1L || !v %in% point_values) {
      fail(sprintf("an unknown %s", f))
    }
    res[[f]] <- v
  }
  n <- length(sf$c)
  if (is.null(res$x)) res$x <- rep(NA_real_, n)
  if (!is.numeric(res$x) || length(res$x) < n) {
    fail(sprintf("`x` of length %d for %d columns", length(res$x), n))
  }
  if (is.null(res$objective)) res$objective <- NA_real_
  res
}
