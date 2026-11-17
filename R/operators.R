# Operator overloading for opt_expr.

new_relation <- function(lhs, rhs, sense) {
  # an infinite numeric side cannot be moved into the expression (y - Inf
  # has no value); it stays a bound, with the expression on the left
  inf_side <- function(v) is.numeric(v) && !is_expr(v) && any(is.infinite(v))
  if (inf_side(lhs)) {
    tmp <- lhs
    lhs <- rhs
    rhs <- tmp
    sense <- switch(sense, "<=" = ">=", ">=" = "<=", "==" = "==")
  }
  if (inf_side(rhs)) {
    return(structure(list(expr = as_expr(lhs), sense = sense,
                          rhs = as.numeric(rhs)), class = "opt_relation"))
  }
  lhs <- as_expr(lhs)
  rhs <- as_expr(rhs)
  e <- expr_add(lhs, rhs, sign = -1)
  structure(list(expr = e, sense = sense), class = "opt_relation")
}

#' Arithmetic and comparison operators for expressions
#'
#' Expressions support `+`, `-`, `*`, `/` and `^` elementwise, with numeric
#' vectors, matrices and other expressions of the same model. Only
#' length-one operands are recycled; other length mismatches are errors.
#' The comparisons `<=`, `>=` and `==` build constraint relations for
#' [opt_constraint()]; they never return `TRUE` or `FALSE`. Strict
#' inequalities `<` and `>` and `!=` are rejected because they do not define
#' closed feasible sets. A power `x^p` whose exponent is an expression of
#' parameters only stays a power node and becomes an ordinary constant power
#' when the parameter values are substituted, so `x^p` with `p = 2` is a
#' square even for negative `x`; an exponent containing decision variables
#' is rewritten as `exp(p * log(x))`, defined for positive `x`. Logical
#' operators are not overloaded; R's `&&` and
#' `||` cannot be overloaded at all, and logical constraints between
#' expressions are not supported.
#'
#' @param e1,e2 Expressions or numeric values.
#' @return An `opt_expr` for arithmetic, an `opt_relation` for comparisons.
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x", dim = 2)
#' 2 * x + 1
#' sum(x^2)
#' x[1] + x[2] <= 4
#' @name opt_operators
NULL

#' @rdname opt_operators
#' @export
Ops.opt_expr <- function(e1, e2) {
  if (missing(e2)) {
    switch(.Generic,
      "-" = return(expr_scale(e1, -1)),
      "+" = return(e1),
      stop(sprintf("Unary `%s` is not defined for expressions.", .Generic),
           call. = FALSE))
  }
  switch(.Generic,
    "+" = expr_add(e1, e2),
    "-" = expr_add(e1, e2, sign = -1),
    "*" = expr_mul(e1, e2),
    "/" = expr_div(e1, e2),
    "^" = expr_pow(e1, e2),
    "<=" = new_relation(e1, e2, "<="),
    ">=" = new_relation(e1, e2, ">="),
    "==" = new_relation(e1, e2, "=="),
    "<" = , ">" = stop(
      "Strict inequalities are not supported; use `<=` or `>=`.",
      call. = FALSE),
    "!=" = stop("`!=` does not define a closed feasible set and is not supported.",
                call. = FALSE),
    "&" = , "|" = , "!" = stop(
      "Logical operators are not defined for optimization expressions.",
      call. = FALSE),
    stop(sprintf("Operator `%s` is not supported for expressions.", .Generic),
         call. = FALSE))
}

#' @rdname opt_operators
#' @param x An expression.
#' @param ... For `log()`, an optional `base`; otherwise unused.
#' @export
Math.opt_expr <- function(x, ...) {
  if (.Generic == "log" && length(list(...))) {
    base <- list(...)[[1L]]
    if (!is.numeric(base) || length(base) != 1L || !is.finite(base) ||
        base <= 0 || base == 1) {
      stop("`base` must be a positive number other than 1.", call. = FALSE)
    }
    return(expr_scale(expr_unary("log", x), 1 / log(base)))
  }
  if (.Generic %in% unary_ops) return(expr_unary(.Generic, x))
  if (.Generic == "cumsum") {
    n <- x$len
    M <- Matrix::sparseMatrix(i = rep.int(seq_len(n), seq_len(n)),
                              j = sequence(seq_len(n)), x = 1, dims = c(n, n))
    return(expr_linmap(M, x))
  }
  if (.Generic %in% c("log2", "log10")) {
    base <- if (.Generic == "log2") 2 else 10
    return(expr_scale(expr_unary("log", x), 1 / log(base)))
  }
  stop(sprintf("`%s()` is not supported for optimization expressions.",
               .Generic), call. = FALSE)
}

#' @rdname opt_operators
#' @param na.rm Ignored; expressions have no missing values.
#' @export
Summary.opt_expr <- function(..., na.rm = FALSE) {
  args <- list(...)
  if (.Generic != "sum") {
    stop(sprintf(paste0("`%s()` of expressions is not supported. ",
                        "Model it with auxiliary variables and linear constraints ",
                        "(an epigraph formulation)."), .Generic), call. = FALSE)
  }
  terms <- lapply(args, function(a) expr_sum(as_expr(a)))
  Reduce(expr_add, terms)
}

#' @rdname opt_operators
#' @export
mean.opt_expr <- function(x, ...) {
  expr_scale(expr_sum(x), 1 / x$len)
}

#' @rdname opt_operators
#' @param y An expression or numeric matrix.
#' @export
matrixOps.opt_expr <- function(x, y) {
  if (.Generic != "%*%") stop("Unsupported matrix operator.", call. = FALSE)
  expr_matmul(x, y)
}

# Matrix shapes of the operands of x %*% y by the rules of base R: a vector
# next to a matrix is a row or a column, whichever conforms; two vectors of
# the same length give their inner product, and a length one vector on the
# left gives an outer product.
matmul_dims <- function(x, y) {
  shape <- function(a) {
    d <- if (is_expr(a)) a$dim else dim(a)
    if (!is.null(d) && length(d) == 2L) as.integer(d) else NULL
  }
  len <- function(a) if (is_expr(a)) a$len else length(a)
  dx <- shape(x)
  dy <- shape(y)
  nx <- len(x)
  ny <- len(y)
  if (is.null(dx) && is.null(dy)) {
    if (nx == ny) return(list(c(1L, nx), c(ny, 1L)))
    if (nx == 1L) return(list(c(1L, 1L), c(1L, ny)))
  } else if (is.null(dx)) {
    if (nx == dy[1L]) return(list(c(1L, nx), dy))
    if (dy[1L] == 1L) return(list(c(nx, 1L), dy))
  } else if (is.null(dy)) {
    if (ny == dx[2L]) return(list(dx, c(ny, 1L)))
    if (dx[2L] == 1L) return(list(dx, c(1L, ny)))
  } else if (dx[2L] == dy[1L]) {
    return(list(dx, dy))
  }
  stop("Non-conformable arguments in %*%.", call. = FALSE)
}

expr_matmul <- function(x, y) {
  if (!is_expr(x) && !is_expr(y)) return(x %*% y)
  dims <- matmul_dims(x, y)
  dx <- dims[[1L]]
  dy <- dims[[2L]]
  if (is_expr(x) && is_expr(y)) {
    if (dx[1L] == 1L && dy[2L] == 1L) {
      # row times column: flatten both so elementwise shapes agree
      return(expr_sum(expr_mul(expr_subset(x, seq_len(x$len)),
                                      expr_subset(y, seq_len(y$len)))))
    }
    # general product: element (r, c) = sum_k x[r, k] * y[k, c]
    r <- rep(seq_len(dx[1L]), times = dy[2L])
    cc <- rep(seq_len(dy[2L]), each = dx[1L])
    nk <- dx[2L]
    xi <- (rep(seq_len(nk), times = length(r)) - 1L) * dx[1L] + rep(r, each = nk)
    yi <- (rep(cc, each = nk) - 1L) * dy[1L] + rep(seq_len(nk), times = length(r))
    prod <- expr_mul(expr_subset(x, xi), expr_subset(y, yi))
    M <- Matrix::sparseMatrix(i = rep(seq_along(r), each = nk),
                              j = seq_along(xi), x = 1,
                              dims = c(length(r), length(xi)))
    return(expr_linmap(M, prod, dim = c(dx[1L], dy[2L])))
  }
  if (is_expr(y)) {
    A <- as_sparse_dims(x, dx)
    M <- if (dy[2L] == 1L) A else
      Matrix::kronecker(Matrix::Diagonal(dy[2L]), A)
    return(expr_linmap(M, y, dim = c(dx[1L], dy[2L])))
  }
  Bt <- Matrix::t(as_sparse_dims(y, dy))
  M <- if (dx[1L] == 1L) Bt else Matrix::kronecker(Bt, Matrix::Diagonal(dx[1L]))
  expr_linmap(M, x, dim = c(dx[1L], dy[2L]))
}

# A numeric vector, matrix or Matrix as a sparse matrix of the given shape.
# A sparse Matrix is never made dense on the way.
as_sparse_dims <- function(a, d) {
  if (methods::is(a, "Matrix")) {
    return(methods::as(methods::as(a, "generalMatrix"), "CsparseMatrix"))
  }
  Matrix::Matrix(matrix(as.vector(a), d[1L], d[2L]), sparse = TRUE,
                 doDiag = FALSE)
}

#' Shape and indexing of expressions
#'
#' Expressions behave like numeric vectors or matrices of scalar
#' expressions. Indexing accepts the same positive, negative, logical and
#' name indices as base R arrays. `c()` concatenates expressions of one
#' model into a vector, and `t()` transposes a matrix expression.
#'
#' @param x An expression.
#' @param i,j Indices.
#' @param drop Drop dimensions of extent one, as in base R.
#' @param ... Further indices or expressions.
#' @param value New dimensions for `dim<-`.
#' @return An expression, or for `length()` and `dim()` integers.
#' @examples
#' m <- opt_model()
#' X <- opt_variable(m, "X", dim = c(2, 3))
#' dim(X)
#' X[1, ]
#' t(X)
#' c(X[1, 1], X[2, 3])
#' @name opt_expr_shape
NULL

#' @rdname opt_expr_shape
#' @export
`[.opt_expr` <- function(x, i, j, ..., drop = TRUE) {
  nidx <- length(sys.call()) - 2L - !missing(drop)
  if (nidx == 1L && is.numeric(i) && !is.matrix(i) && length(i) && !anyNA(i) &&
      all(i >= 1) && all(i <= x$len)) {
    return(expr_subset(x, i))
  }
  d <- x$dim %||% x$len
  dn <- x$dimnames
  if (!is.null(dn) && length(dn) != length(d)) dn <- NULL
  pos <- array(seq_len(x$len), dim = d, dimnames = dn)
  if (nidx <= 1L) {
    idx <- if (missing(i)) pos[] else pos[i]
  } else {
    if (length(d) != 2L) stop("Use one index for vector expressions.",
                              call. = FALSE)
    if (missing(i)) i <- TRUE
    if (missing(j)) j <- TRUE
    idx <- pos[i, j, drop = drop]
  }
  nd <- dim(idx)
  if (!is.null(nd) && length(nd) < 2L) nd <- NULL
  expr_subset(x, as.vector(idx), dim = nd)
}

#' @rdname opt_expr_shape
#' @export
`[[.opt_expr` <- function(x, i, ...) {
  if (!is.character(i)) {
    stop("Use `[` to select elements of an expression.", call. = FALSE)
  }
  .subset2(x, i)
}

#' @rdname opt_expr_shape
#' @export
length.opt_expr <- function(x) x$len

#' @rdname opt_expr_shape
#' @export
dim.opt_expr <- function(x) x$dim

#' @rdname opt_expr_shape
#' @export
`dim<-.opt_expr` <- function(x, value) {
  if (!is.null(value)) {
    value <- as.integer(value)
    if (prod(value) != x$len) stop("dims do not match the length of the expression.",
                                   call. = FALSE)
    if (length(value) < 2L) value <- NULL
  }
  x$dim <- value
  x
}

#' @rdname opt_expr_shape
#' @export
t.opt_expr <- function(x) {
  d <- x$dim %||% c(x$len, 1L)
  if (length(d) != 2L) stop("t() needs a matrix expression.", call. = FALSE)
  idx <- as.vector(t(matrix(seq_len(x$len), d[1L], d[2L])))
  expr_subset(x, idx, dim = rev(d))
}

#' @rdname opt_expr_shape
#' @export
c.opt_expr <- function(...) {
  args <- lapply(list(...), as_expr)
  lens <- vapply(args, `[[`, 0L, "len")
  n <- sum(lens)
  off <- cumsum(c(0L, lens))
  model <- NULL
  for (a in args) {
    if (!is.null(a$model)) {
      if (!is.null(model) && !identical(model, a$model)) {
        stop("Cannot combine expressions that belong to different models.",
             call. = FALSE)
      }
      model <- a$model
    }
  }
  # the polynomial parts are stacked in one pass, with row indices shifted
  # by each argument's offset; nonlinear parts are placed by linear maps
  pl <- vapply(args, function(a) a$type == "poly", NA)
  pa <- args[pl]
  po <- off[-length(off)][pl]
  pull <- function(k, empty) unlist(lapply(pa, .subset2, k),
                                   use.names = FALSE) %||% empty
  nl <- vapply(pa, function(a) length(a$li), 0L)
  nq <- vapply(pa, function(a) length(a$qi), 0L)
  const <- numeric(n)
  for (k in seq_along(pa)) {
    const[po[k] + seq_len(pa[[k]]$len)] <- rep_len(pa[[k]]$const, pa[[k]]$len)
  }
  out <- new_poly(model, n, NULL, const = const,
                  li = pull("li", integer()) + rep.int(po, nl),
                  lj = pull("lj", integer()), lv = pull("lv", numeric()),
                  qi = pull("qi", integer()) + rep.int(po, nq),
                  qj = pull("qj", integer()), qk = pull("qk", integer()),
                  qv = pull("qv", numeric()))
  for (k in which(!pl)) {
    a <- args[[k]]
    M <- Matrix::sparseMatrix(i = off[k] + seq_len(a$len),
                              j = seq_len(a$len), x = 1, dims = c(n, a$len))
    out <- expr_add(out, expr_linmap(M, a))
  }
  out
}

#' @export
rep.opt_expr <- function(x, ...) {
  idx <- rep(seq_len(x$len), ...)
  expr_subset(x, idx)
}

# S4 methods so that Matrix objects on the left dispatch to optcompR.
#' @importFrom methods setOldClass setMethod signature
NULL

methods::setOldClass("opt_expr")

#' @rdname opt_operators
#' @export
setMethod("%*%", signature("Matrix", "opt_expr"), function(x, y) {
  expr_matmul(x, y)
})

#' @rdname opt_operators
#' @export
setMethod("%*%", signature("opt_expr", "Matrix"), function(x, y) {
  expr_matmul(x, y)
})
