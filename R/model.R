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
