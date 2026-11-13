valid_types <- c("continuous", "integer", "binary")

normalize_dim <- function(dim, dimnames) {
  if (!is.null(dimnames)) {
    if (!is.list(dimnames)) dimnames <- list(dimnames)
    dn_dim <- vapply(dimnames, length, 0L)
    if (missing(dim) || is.null(dim)) dim <- dn_dim
    if (length(dim) != length(dimnames) || any(dn_dim != dim)) {
      stop("`dimnames` does not match `dim`.", call. = FALSE)
    }
    dimnames <- lapply(dimnames, as.character)
  }
  if (!is.numeric(dim) || anyNA(dim) || any(dim < 1) ||
      any(dim != round(dim))) {
    stop("`dim` must be a vector of positive integers.", call. = FALSE)
  }
  dim <- as.integer(dim)
  list(dim = dim, dimnames = dimnames, len = prod(dim),
       expr_dim = if (length(dim) > 1L) dim else NULL)
}

register_block <- function(model, name, kind, labels, dims) {
  check_scalar_string(name, "name")
  if (!is.na(index_get(model$block_index, name))) {
    stop(sprintf("The model already has a variable or parameter named '%s'.",
                 name), call. = FALSE)
  }
  clash <- labels[labels %in% model$sym_name | duplicated(labels)]
  if (length(clash)) {
    stop(sprintf("The element name '%s' is already used in this model.",
                 clash[1L]), call. = FALSE)
  }
  n0 <- n_symbols(model)
  sym <- n0 + seq_len(dims$len)
  bid <- length(model$blocks) + 1L
  model$blocks[[bid]] <- list(name = name, kind = kind, sym = sym,
                              dim = dims$dim, dimnames = dims$dimnames)
  model$block_index[[name]] <- bid
  model$sym_kind <- c(model$sym_kind, rep.int(kind, dims$len))
  model$sym_name <- c(model$sym_name, labels)
  model$sym_block <- c(model$sym_block, rep.int(bid, dims$len))
  sym
}

symbol_expr <- function(model, sym, dim) {
  n <- length(sym)
  new_poly(model, n, dim, const = numeric(n), li = seq_len(n), lj = sym,
           lv = rep.int(1, n))
}

check_bounds <- function(lower, upper, type, labels) {
  check_numeric(lower, "lower")
  check_numeric(upper, "upper")
  bad <- lower > upper
  if (any(bad)) {
    stop(sprintf("Contradictory bounds for %s: lower %s > upper %s.",
                 labels[which(bad)[1L]], fmt_num(lower[bad][1L]),
                 fmt_num(upper[bad][1L])), call. = FALSE)
  }
  if (any(lower == Inf) || any(upper == -Inf)) {
    stop("A lower bound of Inf or an upper bound of -Inf leaves no feasible value.",
         call. = FALSE)
  }
  bin <- type == "binary"
  if (any(bin & (lower > 1 | upper < 0))) {
    stop("Binary variables need bounds that intersect [0, 1].", call. = FALSE)
  }
  intlike <- type %in% c("integer", "binary")
  if (any(intlike & is.finite(lower) & is.finite(upper) &
          ceiling(lower - 1e-9) > floor(upper + 1e-9))) {
    stop("An integer variable has bounds that contain no integer.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Add decision variables to a model
#'
#' Creates a scalar, vector, or array block of decision variables and returns
#' an expression that refers to them. Variables can be continuous, integer or
#' binary. Bounds are validated: contradictory bounds, integer bounds without
#' an integer between them, and binary bounds outside `[0, 1]` are errors.
#'
#' @param model An [opt_model()].
#' @param name Block name. Must be unique among the model's variables and
#'   parameters. Element names are `name[i]` or `name[i,j]`, or use
#'   `dimnames`.
#' @param dim Integer dimensions: `1` for a scalar, `n` for a vector,
#'   `c(n, p)` for a matrix. Taken from `dimnames` when omitted.
#' @param type `"continuous"`, `"integer"` or `"binary"`.
#' @param lower,upper Bounds, recycled from length one. Binary variables have
#'   bounds within `[0, 1]`; the defaults for other types are free.
#' @param dimnames Optional character labels, a vector for one dimension or a
#'   list with one entry per dimension.
#' @param start Optional start values (used by solvers that accept warm
#'   starts).
#' @return An `opt_expr` with the shape given by `dim`.
#' @seealso [opt_set_bounds()], [opt_set_type()], [opt_fix()]
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x", lower = 0)
#' y <- opt_variable(m, "y", dim = 3, type = "integer", lower = 0, upper = 10)
#' z <- opt_variable(m, "z", type = "binary",
#'                   dimnames = list(c("a", "b"), c("u", "v")))
#' opt_variables(m)
#' @export
opt_variable <- function(model, name, dim = 1L,
                         type = c("continuous", "integer", "binary"),
                         lower = NULL, upper = NULL, dimnames = NULL,
                         start = NULL) {
  check_model(model)
  type <- match.arg(type)
  dims <- normalize_dim(if (missing(dim) && !is.null(dimnames)) NULL else dim,
                        dimnames)
  n <- dims$len
  if (is.null(lower)) lower <- if (type == "binary") 0 else -Inf
  if (is.null(upper)) upper <- if (type == "binary") 1 else Inf
  lower <- recycle_to(as.numeric(lower), n, "lower")
  upper <- recycle_to(as.numeric(upper), n, "upper")
  labels <- element_labels(name, n, if (length(dims$dim) > 1L || n > 1L)
    dims$dim else NULL, dims$dimnames)
  types <- rep.int(type, n)
  check_bounds(lower, upper, types, labels)
  if (type == "binary") {
    lower <- pmax(lower, 0)
    upper <- pmin(upper, 1)
  }
  st <- if (is.null(start)) rep(NA_real_, n) else
    recycle_to(as.numeric(start), n, "start")
  sym <- register_block(model, name, "var", labels, dims)
  model$lower <- c(model$lower, lower)
  model$upper <- c(model$upper, upper)
  model$type <- c(model$type, types)
  model$start <- c(model$start, st)
  model$value <- c(model$value, rep(NA_real_, n))
  touch_structure(model)
  out <- symbol_expr(model, sym, dims$expr_dim)
  out$dimnames <- dims$dimnames
  out
}

#' Add parameters to a model
#'
#' Parameters are fixed numbers that enter expressions symbolically, so the
#' model can be solved again for new data with [opt_set_parameter()] without
#' rebuilding its constraints. Problem classification and convexity checks
#' happen after the current values are substituted, so a parameter that
#' multiplies a quadratic term still gives a quadratic program.
#'
#' @inheritParams opt_variable
#' @param value Numeric values, recycled from length one.
#' @param dim Dimensions; by default the length (or `dim()`) of `value`.
#' @return An `opt_expr` referring to the parameters.
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x")
#' p <- opt_parameter(m, "p", value = 1)
#' opt_objective(m, (x - p)^2)
#' opt_value(opt_solve(m), x)
#' opt_set_parameter(m, p, 5)
#' opt_value(opt_solve(m), x)
#' @export
opt_parameter <- function(model, name, value, dim = NULL, dimnames = NULL) {
  check_model(model)
  if (is.null(dim) && is.null(dimnames)) dim <- if (is.null(base::dim(value)))
    length(value) else base::dim(value)
  dims <- normalize_dim(dim, dimnames)
  n <- dims$len
  value <- recycle_to(as.numeric(value), n, "value")
  check_numeric(value, "value")
  if (any(!is.finite(value))) stop("Parameter values must be finite.",
                                   call. = FALSE)
  labels <- element_labels(name, n, if (n > 1L) dims$dim else NULL,
                           dims$dimnames)
  sym <- register_block(model, name, "par", labels, dims)
  model$lower <- c(model$lower, rep(NA_real_, n))
  model$upper <- c(model$upper, rep(NA_real_, n))
  model$type <- c(model$type, rep(NA_character_, n))
  model$start <- c(model$start, rep(NA_real_, n))
  model$value <- c(model$value, value)
  touch_structure(model)
  touch_parameters(model)
  out <- symbol_expr(model, sym, dims$expr_dim)
  out$dimnames <- dims$dimnames
  out
}

# Resolve an expression (or block name) to the symbol indices it denotes.
# Only plain symbol expressions, such as variables or slices of them, are
# accepted.
resolve_symbols <- function(model, x, kind) {
  if (is.character(x)) {
    check_scalar_string(x, "x")
    bid <- index_get(model$block_index, x)
    if (is.na(bid)) stop(sprintf("No variable or parameter named '%s'.", x),
                         call. = FALSE)
    sym <- model$blocks[[bid]]$sym
  } else {
    if (!is_expr(x) || x$type != "poly" || !identical(x$model, model)) {
      stop("`x` must be a variable or parameter of this model, or its name.",
           call. = FALSE)
    }
    x <- compact_poly(x)
    ok <- length(x$qv) == 0L && all(x$const == 0) &&
      length(x$lv) == x$len && all(x$lv == 1) &&
      identical(sort(x$li), seq_len(x$len))
    if (!ok) stop("`x` must be a plain variable or parameter, not an expression.",
                  call. = FALSE)
    sym <- x$lj[order(x$li)]
  }
  if (!all(model$sym_kind[sym] == kind)) {
    stop(sprintf("`x` must refer to %s.",
                 if (kind == "var") "decision variables" else "parameters"),
         call. = FALSE)
  }
  sym
}

#' Change parameter values
#'
#' Updates parameter values in place. Cached standard forms built for the old
#' values are invalidated; the next [opt_solve()] rebuilds the numerical
#' standard form from the unchanged expressions (the symbolic model is kept;
#' the lowering and the backend model are built again).
#'
#' @param model An [opt_model()].
#' @param x A parameter expression returned by [opt_parameter()] (or a slice
#'   of it), or the parameter block name.
#' @param value New values, recycled from length one.
#' @return The model, invisibly.
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x", lower = 0)
#' b <- opt_parameter(m, "b", value = 2)
#' opt_constraint(m, x >= b)
#' opt_objective(m, x)
#' opt_objective_value(opt_solve(m))
#' opt_set_parameter(m, b, 7)
#' opt_objective_value(opt_solve(m))
#' @export
opt_set_parameter <- function(model, x, value) {
  check_model(model)
  sym <- resolve_symbols(model, x, "par")
  value <- recycle_to(as.numeric(value), length(sym), "value")
  check_numeric(value, "value")
  if (any(!is.finite(value))) stop("Parameter values must be finite.",
                                   call. = FALSE)
  model$value[sym] <- value
  touch_parameters(model)
  invisible(model)
}

#' Change variable bounds
#'
#' @param model An [opt_model()].
#' @param x A variable expression (or slice) or the variable block name.
#' @param lower,upper New bounds, recycled from length one; `NULL` keeps the
#'   current bound.
#' @return The model, invisibly.
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x", dim = 2)
#' opt_set_bounds(m, x, lower = 0, upper = c(1, 2))
#' opt_variables(m)
#' @export
opt_set_bounds <- function(model, x, lower = NULL, upper = NULL) {
  check_model(model)
  sym <- resolve_symbols(model, x, "var")
  n <- length(sym)
  lo <- if (is.null(lower)) model$lower[sym] else
    recycle_to(as.numeric(lower), n, "lower")
  up <- if (is.null(upper)) model$upper[sym] else
    recycle_to(as.numeric(upper), n, "upper")
  check_bounds(lo, up, model$type[sym], model$sym_name[sym])
  model$lower[sym] <- lo
  model$upper[sym] <- up
  touch_parameters(model)
  invisible(model)
}

#' Fix variables to values
#'
#' Sets both bounds to `value`.
#'
#' @inheritParams opt_set_bounds
#' @param value Values, recycled from length one.
#' @return The model, invisibly.
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x", lower = 0)
#' opt_fix(m, x, 3)
#' opt_variables(m)
#' @export
opt_fix <- function(model, x, value) {
  opt_set_bounds(model, x, lower = value, upper = value)
}

#' Change variable types
#'
#' Changing a type changes the problem class, for example from a linear to a
#' mixed-integer linear program. Setting `"binary"` also restricts the bounds
#' to `[0, 1]`.
#'
#' @inheritParams opt_set_bounds
#' @param type `"continuous"`, `"integer"` or `"binary"`.
#' @return The model, invisibly.
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x", lower = 0, upper = 5)
#' opt_set_type(m, x, "integer")
#' opt_variables(m)
#' @export
opt_set_type <- function(model, x, type = c("continuous", "integer", "binary")) {
  check_model(model)
  type <- match.arg(type)
  sym <- resolve_symbols(model, x, "var")
  lo <- model$lower[sym]
  up <- model$upper[sym]
  if (type == "binary") {
    lo <- pmax(lo, 0)
    up <- pmin(up, 1)
  }
  check_bounds(lo, up, rep.int(type, length(sym)), model$sym_name[sym])
  model$type[sym] <- type
  model$lower[sym] <- lo
  model$upper[sym] <- up
  touch_structure(model)
  invisible(model)
}

#' Set start values for variables
#'
#' Start values are passed to solvers that accept a warm start (for example
#' HiGHS for mixed-integer problems and NLopt for nonlinear problems).
#'
#' @inheritParams opt_set_bounds
#' @param value Start values, recycled from length one; `NA` clears them.
#' @return The model, invisibly.
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x", dim = 2)
#' opt_set_start(m, x, c(1, 2))
#' @export
opt_set_start <- function(model, x, value) {
  check_model(model)
  sym <- resolve_symbols(model, x, "var")
  model$start[sym] <- recycle_to(as.numeric(value), length(sym), "value")
  touch_parameters(model)
  invisible(model)
}

#' Retrieve a variable or parameter block by name
#'
#' @param model An [opt_model()].
#' @param name Block name.
#' @return The expression for the whole block.
#' @examples
#' m <- opt_model()
#' opt_variable(m, "x", dim = 3)
#' opt_get(m, "x")
#' @export
opt_get <- function(model, name) {
  check_model(model)
  check_scalar_string(name, "name")
  bid <- index_get(model$block_index, name)
  if (is.na(bid)) stop(sprintf("No variable or parameter named '%s'.", name),
                       call. = FALSE)
  b <- model$blocks[[bid]]
  out <- symbol_expr(model, b$sym, if (length(b$dim) > 1L) b$dim else NULL)
  out$dimnames <- b$dimnames
  out
}
