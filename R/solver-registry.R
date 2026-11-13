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
  if (cls$class %in% e$convex_only && cls$convexity == "nonconvex") {
    return(sprintf("%s requires a convex %s, and this model is nonconvex",
                   e$name, cls$class))
  }
  NULL
}

select_solver <- function(cls, solver) {
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
    return(e)
  }
  stop(sprintf("No installed solver supports this %s (%s). %s.", cls$class,
               cls$convexity, paste(unique(reasons), collapse = "; ")),
       call. = FALSE)
}
