#' Add constraints to a model
#'
#' Adds a block of constraints. The block is either a relation built with
#' `<=`, `>=` or `==` (for example `A %*% x <= b`), or an expression with
#' explicit `lower` and `upper` bounds, which gives ranged rows
#' `lower <= expr <= upper`. Vector relations add one row per element.
#'
#' Constraints are stored exactly as written. They are never approximated;
#' a nonlinear constraint makes the model nonlinear, and solvers that cannot
#' handle it are not selected.
#'
#' @param model An [opt_model()].
#' @param relation An `opt_relation` from a comparison, or an `opt_expr`
#'   when `lower` or `upper` is given.
#' @param name Block name, unique within the model. Defaults to `c1`, `c2`,
#'   and so on.
#' @param lower,upper Row bounds for the ranged form, recycled from length
#'   one.
#' @return An `opt_constraint` reference, used by [opt_dual()],
#'   [opt_slack()] and [opt_delete_constraint()]. Returned invisibly.
#' @examples
#' m <- opt_model("max")
#' x <- opt_variable(m, "x", dim = 2, lower = 0)
#' cap <- opt_constraint(m, x[1] + x[2] <= 4, name = "capacity")
#' opt_constraint(m, x[1] - x[2], lower = -1, upper = 1, name = "balance")
#' opt_objective(m, 3 * x[1] + 2 * x[2])
#' r <- opt_solve(m)
#' opt_dual(r, cap)
#' @export
opt_constraint <- function(model, relation, name = NULL, lower = NULL,
                           upper = NULL) {
  check_model(model)
  if (inherits(relation, "opt_relation")) {
    if (!is.null(lower) || !is.null(upper)) {
      stop("Give either a relation or an expression with `lower`/`upper`, not both.",
           call. = FALSE)
    }
    e <- relation$expr
    n <- e$len
    lo <- switch(relation$sense, "<=" = rep(-Inf, n), ">=" = rep(0, n),
                 "==" = rep(0, n))
    up <- switch(relation$sense, "<=" = rep(0, n), ">=" = rep(Inf, n),
                 "==" = rep(0, n))
    sense <- relation$sense
  } else if (is_expr(relation)) {
    if (is.null(lower) && is.null(upper)) {
      stop("An expression constraint needs `lower` and/or `upper`.",
           call. = FALSE)
    }
    e <- relation
    n <- e$len
    lo <- recycle_to(as.numeric(lower %||% -Inf), n, "lower")
    up <- recycle_to(as.numeric(upper %||% Inf), n, "upper")
    if (anyNA(lo) || anyNA(up) || any(lo > up)) {
      stop("Constraint bounds must satisfy lower <= upper.", call. = FALSE)
    }
    sense <- "range"
  } else if (is.logical(relation)) {
    stop("The constraint evaluated to TRUE/FALSE; it contains no optimization expression.",
         call. = FALSE)
  } else {
    stop("`relation` must be a comparison such as `x + y <= 4`.",
         call. = FALSE)
  }
  if (is.null(e$model)) {
    stop("The constraint contains no variables or parameters.", call. = FALSE)
  }
  if (!identical(e$model, model)) {
    stop("The constraint belongs to a different model.", call. = FALSE)
  }
  if (is.null(name)) {
    k <- length(model$cons) + 1L
    repeat {
      name <- paste0("c", k)
      if (is.na(index_get(model$cons_index, name))) break
      k <- k + 1L
    }
  }
  check_scalar_string(name, "name")
  if (!is.na(index_get(model$cons_index, name))) {
    stop(sprintf("The model already has a constraint block named '%s'.", name),
         call. = FALSE)
  }
  bid <- length(model$cons) + 1L
  cons <- model$cons
  model$cons <- NULL # drop the reference so the list grows in place
  cons[[bid]] <- list(name = name, expr = e, lower = lo, upper = up,
                      sense = sense)
  model$cons <- cons
  model$cons_index[[name]] <- bid
  touch_structure(model)
  invisible(structure(list(model = model, name = name),
                      class = "opt_constraint"))
}

#' @export
print.opt_constraint <- function(x, ...) {
  b <- constraint_block(x$model, x$name)
  cat(sprintf("<opt_constraint '%s': %d row%s, %s>\n", x$name, b$expr$len,
              if (b$expr$len == 1L) "" else "s", expr_kind_label(b$expr)))
  invisible(x)
}

constraint_block <- function(model, name) {
  bid <- index_get(model$cons_index, name)
  if (is.na(bid)) stop(sprintf("No constraint block named '%s'.", name),
                       call. = FALSE)
  model$cons[[bid]]
}

#' Delete a constraint block
#'
#' @param model An [opt_model()].
#' @param constraint An `opt_constraint` reference or the block name.
#' @return The model, invisibly.
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x", lower = 0)
#' cn <- opt_constraint(m, x >= 1)
#' opt_delete_constraint(m, cn)
#' opt_constraints(m)
#' @export
opt_delete_constraint <- function(model, constraint) {
  check_model(model)
  name <- if (inherits(constraint, "opt_constraint")) constraint$name else
    constraint
  check_scalar_string(name, "constraint")
  bid <- index_get(model$cons_index, name)
  if (is.na(bid)) stop(sprintf("No constraint block named '%s'.", name),
                       call. = FALSE)
  model$cons[[bid]] <- NULL
  nm <- vapply(model$cons, `[[`, "", "name")
  model$cons_index <- list2env(as.list(stats::setNames(seq_along(nm), nm)),
                               parent = emptyenv())
  touch_structure(model)
  invisible(model)
}

#' Set the objective of a model
#'
#' @param model An [opt_model()].
#' @param expr A scalar expression. Constant terms are kept and reported in
#'   objective values.
#' @param sense Optional new sense, `"min"` or `"max"`; by default the
#'   model's sense is kept.
#' @return The model, invisibly.
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x", lower = 1)
#' opt_objective(m, 2 * x + 5)
#' opt_objective_value(opt_solve(m))
#' @export
opt_objective <- function(model, expr, sense = NULL) {
  check_model(model)
  expr <- as_expr(expr)
  if (expr$len != 1L) {
    stop(sprintf("The objective must be a scalar expression, not length %d. Use sum().",
                 expr$len), call. = FALSE)
  }
  if (!is.null(expr$model) && !identical(expr$model, model)) {
    stop("The objective belongs to a different model.", call. = FALSE)
  }
  expr$model <- model
  if (!is.null(sense)) model$sense <- match.arg(sense, c("min", "max"))
  model$objective <- expr
  touch_structure(model)
  invisible(model)
}

#' @export
print.opt_relation <- function(x, ...) {
  cat(sprintf("<opt_relation: %d row%s, %s %s 0>\n", x$expr$len,
              if (x$expr$len == 1L) "" else "s", expr_kind_label(x$expr),
              x$sense))
  e <- x$expr
  n <- min(e$len, 10L)
  if (e$type == "poly") {
    e <- compact_poly(e)
    lhs <- vapply(seq_len(n), format_poly_element, "", e = e,
                  drop_const = TRUE)
    rhs <- fmt_num(-e$const[seq_len(n)])
  } else {
    lhs <- format_expr(e, n)
    rhs <- "0"
  }
  cat(paste0("  ", lhs, " ", x$sense, " ", rhs), sep = "\n")
  if (e$len > n) cat(sprintf("  ... %d more\n", e$len - n))
  invisible(x)
}
