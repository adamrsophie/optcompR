.optcompR_state <- new.env(parent = emptyenv())
.optcompR_state$counter <- 0L

new_model_id <- function() {
  .optcompR_state$counter <- .optcompR_state$counter + 1L
  # a counter only: drawing random numbers would change the user's RNG
  # stream; identity is the environment itself
  sprintf("m%d", .optcompR_state$counter)
}

#' Create an optimization model
#'
#' An `opt_model` is a mutable container for decision variables, parameters,
#' constraints and one objective. Models have reference semantics: functions
#' such as [opt_variable()] and [opt_constraint()] change the model in place,
#' and every copy of the object refers to the same model. Expressions remember
#' the model they were built from, and combining expressions from two models is
#' an error.
#'
#' @section Lifetime and serialization:
#' A model lives as long as some object refers to it. Saving a model with
#' [saveRDS()] stores its variables, constraints and objective. Expressions
#' saved in a separate file refer to their own copy of the model after
#' loading, so save a model together with the expressions you want to reuse,
#' for example in one list.
#'
#' @param sense Optimization sense, `"min"` or `"max"`.
#' @param name Optional model name used in printing and reports.
#' @return An object of class `opt_model`.
#' @seealso [opt_variable()], [opt_constraint()], [opt_objective()],
#'   [opt_solve()]
#' @examples
#' m <- opt_model("max")
#' x <- opt_variable(m, "x", lower = 0, upper = 2)
#' y <- opt_variable(m, "y", lower = 0, upper = 3)
#' opt_constraint(m, x + y <= 4, name = "capacity")
#' opt_objective(m, 3 * x + 2 * y)
#' m
#' @export
opt_model <- function(sense = c("min", "max"), name = NULL) {
  sense <- match.arg(sense)
  if (!is.null(name)) check_scalar_string(name, "name")
  m <- new.env(parent = emptyenv())
  m$id <- new_model_id()
  m$name <- name
  m$sense <- sense
  m$sym_kind <- character()
  m$sym_name <- character()
  m$sym_block <- integer()
  m$lower <- numeric()
  m$upper <- numeric()
  m$type <- character()
  m$start <- numeric()
  m$value <- numeric()
  m$blocks <- list()
  m$block_index <- new.env(parent = emptyenv())
  m$cons <- list()
  m$cons_index <- new.env(parent = emptyenv())
  m$objective <- NULL
  m$structure_version <- 0L
  m$parameter_version <- 0L
  m$cache <- new.env(parent = emptyenv())
  class(m) <- "opt_model"
  m
}

check_model <- function(model) {
  if (!inherits(model, "opt_model")) {
    stop("`model` must be an `opt_model` created by opt_model().",
         call. = FALSE)
  }
  invisible(model)
}

touch_structure <- function(model) {
  model$structure_version <- model$structure_version + 1L
  invisible(model)
}

# Evaluate `code`, which adds several variables and constraints to `model`;
# if it fails, the model is put back as it was, so a rejected call leaves no
# partial construction behind. The versions move forward (and the caches are
# cleared) so that no cached lowering can be mistaken for the restored state.
atomically <- function(model, code) {
  snap <- as.list.environment(model, all.names = TRUE)
  bi <- as.list(model$block_index, all.names = TRUE)
  ci <- as.list(model$cons_index, all.names = TRUE)
  ok <- FALSE
  on.exit(if (!ok) {
    sv <- model$structure_version
    pv <- model$parameter_version
    rm(list = setdiff(ls(model, all.names = TRUE), names(snap)), envir = model)
    for (k in setdiff(names(snap), c("block_index", "cons_index", "cache"))) {
      assign(k, snap[[k]], envir = model)
    }
    rm(list = ls(model$block_index, all.names = TRUE), envir = model$block_index)
    rm(list = ls(model$cons_index, all.names = TRUE), envir = model$cons_index)
    list2env(bi, model$block_index)
    list2env(ci, model$cons_index)
    model$cache <- new.env(parent = emptyenv())
    model$structure_version <- max(sv, snap$structure_version) + 1L
    model$parameter_version <- max(pv, snap$parameter_version) + 1L
  }, add = TRUE)
  out <- force(code)
  ok <- TRUE
  out
}

touch_parameters <- function(model) {
  model$parameter_version <- model$parameter_version + 1L
  invisible(model)
}

# Name lookups in O(1): NA when absent.
index_get <- function(env, name) {
  v <- env[[name]]
  if (is.null(v)) NA_integer_ else v
}

n_symbols <- function(model) length(model$sym_kind)

var_symbols <- function(model) which(model$sym_kind == "var")

#' @export
print.opt_model <- function(x, ...) {
  nv <- sum(x$sym_kind == "var")
  np <- sum(x$sym_kind == "par")
  nint <- sum(x$type %in% c("integer", "binary"), na.rm = TRUE)
  nrow <- sum(vapply(x$cons, function(b) b$expr$len, 0L))
  cat("<opt_model", if (!is.null(x$name)) paste0(" '", x$name, "'"), ">\n",
      sep = "")
  cat(sprintf("  sense:       %s\n", x$sense))
  cat(sprintf("  variables:   %d (%d integer or binary)\n", nv, nint))
  if (np > 0L) cat(sprintf("  parameters:  %d\n", np))
  cat(sprintf("  constraints: %d rows in %d blocks\n", nrow, length(x$cons)))
  cat(sprintf("  objective:   %s\n",
              if (is.null(x$objective)) "none (feasibility problem)"
              else expr_kind_label(x$objective)))
  invisible(x)
}

#' List the variables of a model
#'
#' @param model An [opt_model()].
#' @return A data frame with one row per scalar variable: name, block, type,
#'   lower and upper bound, and start value.
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x", dim = 2, lower = 0)
#' opt_variables(m)
#' @export
opt_variables <- function(model) {
  check_model(model)
  s <- var_symbols(model)
  data.frame(
    name = model$sym_name[s],
    block = vapply(model$blocks[model$sym_block[s]], `[[`, "", "name"),
    type = model$type[s],
    lower = model$lower[s],
    upper = model$upper[s],
    start = model$start[s],
    stringsAsFactors = FALSE
  )
}

#' List the constraint blocks of a model
#'
#' @param model An [opt_model()].
#' @return A data frame with one row per constraint block: name, number of
#'   rows, and the polynomial degree or nonlinear kind of its expression.
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x", dim = 3, lower = 0)
#' opt_constraint(m, sum(x) <= 1, name = "budget")
#' opt_constraints(m)
#' @export
opt_constraints <- function(model) {
  check_model(model)
  if (length(model$cons) == 0L) {
    return(data.frame(name = character(), rows = integer(),
                      kind = character(), stringsAsFactors = FALSE))
  }
  data.frame(
    name = vapply(model$cons, `[[`, "", "name"),
    rows = vapply(model$cons, function(b) b$expr$len, 0L),
    kind = vapply(model$cons, function(b) expr_kind_label(b$expr), ""),
    stringsAsFactors = FALSE
  )
}

#' Copy a model
#'
#' Models have reference semantics, so assigning a model to a new name does
#' not copy it. `opt_copy()` returns an independent model with the same
#' variables, parameters, constraints, special constraints, objective and
#' values; changing one model leaves the other unchanged. Expressions built
#' from the original still refer to the original, so retrieve the blocks of
#' the copy with [opt_get()].
#'
#' @param model An [opt_model()].
#' @param name Optional name of the copy; by default the original name.
#' @return A new `opt_model`.
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x", lower = 0)
#' opt_objective(m, x)
#' m2 <- opt_copy(m)
#' opt_set_bounds(m2, opt_get(m2, "x"), lower = 2)
#' opt_objective_value(opt_solve(m))   # 0
#' opt_objective_value(opt_solve(m2))  # 2
#' @export
opt_copy <- function(model, name = NULL) {
  check_model(model)
  if (!is.null(name)) check_scalar_string(name, "name")
  out <- new.env(parent = emptyenv())
  for (k in ls(model, all.names = TRUE)) assign(k, get(k, envir = model), out)
  out$id <- new_model_id()
  if (!is.null(name)) out$name <- name
  out$block_index <- list2env(as.list(model$block_index), parent = emptyenv())
  out$cons_index <- list2env(as.list(model$cons_index), parent = emptyenv())
  out$cache <- new.env(parent = emptyenv())
  class(out) <- "opt_model"
  # expressions stored in the model point at their model: repoint them
  for (k in c("blocks", "cons", "objective", "special")) {
    if (!is.null(out[[k]])) assign(k, relink(out[[k]], model, out), out)
  }
  out
}

# Replace references to the model environment `from` by `to` in x.
relink <- function(x, from, to) {
  if (is.environment(x)) return(if (identical(x, from)) to else x)
  if (is.list(x)) {
    at <- attributes(x)
    x <- lapply(unclass(x), relink, from, to)
    attributes(x) <- at
  }
  x
}

# A new model defining the problem a result was computed for: a copy of its
# model with the constraints, objective, special constraints, sense,
# parameter values, bounds and types the result recorded. The live model is
# not touched.
model_from_result <- function(result) {
  live <- result$model
  m <- opt_copy(live)
  m$cons <- relink(result$cons, live, m)
  nm <- vapply(m$cons, `[[`, "", "name")
  m$cons_index <- list2env(as.list(stats::setNames(seq_along(nm), nm)),
                           parent = emptyenv())
  m$objective <- relink(result$objective, live, m)
  m$special <- relink(result$special, live, m)
  m$sense <- result$sense
  sf <- result$standard_form
  ns <- length(result$z)
  par <- which(m$sym_kind[seq_len(ns)] == "par")
  m$value[par] <- result$z[par]
  vs <- sf$var_sym
  n0 <- length(vs)
  m$lower[vs] <- sf$lower[seq_len(n0)]
  m$upper[vs] <- sf$upper[seq_len(n0)]
  m$type[vs] <- sf$type[seq_len(n0)]
  touch_structure(m)
  touch_parameters(m)
  m
}
