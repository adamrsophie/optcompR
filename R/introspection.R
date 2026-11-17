# Public views of the internal expression representation.

#' Inspect the terms of an expression
#'
#' `opt_terms()` returns the polynomial structure of an expression:
#' its constant, linear and quadratic coefficients per element, using the
#' names of the model's variables and parameters. `opt_degree()` returns the
#' polynomial degree of each element (0, 1 or 2, and `NA` for elements with
#' nonlinear terms; parameters count like variables), and `opt_symbols()`
#' the names of the variables and
#' parameters an expression depends on.
#'
#' @param x An `opt_expr` or a number.
#' @return For `opt_terms()`, a list with `constant` (one value per
#'   element), `linear` (a data frame with columns `element`, `symbol` and
#'   `coef`), `quadratic` (columns `element`, `symbol1`, `symbol2`, `coef`,
#'   meaning `coef * symbol1 * symbol2`) and `nonlinear` (`TRUE` when the
#'   expression has nonlinear nodes, which `opt_terms()` cannot list). For
#'   `opt_degree()` an integer vector; for `opt_symbols()` a character
#'   vector.
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x", dim = 2)
#' e <- c(3 * x[1] + 2, x[1] * x[2] - x[2])
#' opt_terms(e)
#' opt_degree(c(e, exp(x[1])))
#' opt_symbols(e[2])
#' @export
opt_terms <- function(x) {
  x <- as_expr(x)
  if (x$type != "poly") {
    return(list(constant = rep(NA_real_, x$len),
                linear = data.frame(element = integer(), symbol = character(),
                                    coef = numeric()),
                quadratic = data.frame(element = integer(),
                                       symbol1 = character(),
                                       symbol2 = character(),
                                       coef = numeric()),
                nonlinear = TRUE))
  }
  e <- compact_poly(x)
  nm <- function(j) if (is.null(e$model)) character() else e$model$sym_name[j]
  list(constant = e$const,
       linear = data.frame(element = as.integer(e$li), symbol = nm(e$lj),
                           coef = e$lv, stringsAsFactors = FALSE),
       quadratic = data.frame(element = as.integer(e$qi), symbol1 = nm(e$qj),
                              symbol2 = nm(e$qk), coef = e$qv,
                              stringsAsFactors = FALSE),
       nonlinear = FALSE)
}

#' @rdname opt_terms
#' @export
opt_degree <- function(x) {
  x <- as_expr(x)
  vapply(seq_len(x$len), function(k) elem_degree(x, k), 0L)
}

# Polynomial degree of element k, NA when it is not a polynomial (a
# linear map keeps degrees, products add them).
elem_degree <- function(e, k) {
  if (e$type == "poly") {
    d <- 0L
    if (any(e$li == k & e$lv != 0)) d <- 1L
    if (any(e$qi == k & e$qv != 0)) d <- 2L
    return(d)
  }
  arg <- function(a, kk = k) elem_degree(a, if (a$len == 1L) 1L else kk)
  switch(e$op,
    linmap = {
      row <- e$data$M[k, ]
      nz <- which(row != 0)
      if (!length(nz)) 0L else max(vapply(nz, function(t) arg(e$args[[1L]], t), 0L))
    },
    "+" = max(arg(e$args[[1L]]), arg(e$args[[2L]])),
    "*" = arg(e$args[[1L]]) + arg(e$args[[2L]]),
    "^" = {
      p <- rep_len(e$data$p, e$len)[k]
      d <- arg(e$args[[1L]])
      if (identical(d, 0L)) 0L else
        if (!is.na(d) && p >= 0 && p == round(p)) as.integer(d * p) else NA_integer_
    },
    {
      ds <- vapply(e$args, arg, 0L)
      if (all(ds == 0L, na.rm = FALSE) && !anyNA(ds)) 0L else NA_integer_
    })
}

#' @rdname opt_terms
#' @export
opt_symbols <- function(x) {
  x <- as_expr(x)
  if (is.null(x$model)) return(character())
  x$model$sym_name[expr_symbols(x)]
}

#' Evaluate an expression at given values
#'
#' @param x An `opt_expr`.
#' @param values Values of the model's variables: a numeric vector with one
#'   value per variable (in the order of [opt_variables()]) or a named vector
#'   naming every variable. Parameters take their current values.
#' @return A numeric vector (or matrix, for matrix expressions).
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x", dim = 2)
#' opt_evaluate(x[1]^2 + 3 * x[2], c(2, 1))
#' @export
opt_evaluate <- function(x, values) {
  x <- as_expr(x)
  if (is.null(x$model)) {
    out <- x$const
  } else {
    p <- point_vector(x$model, values)
    out <- eval_expr(x, p$z)
  }
  if (!is.null(x$dim)) dim(out) <- x$dim
  out
}
