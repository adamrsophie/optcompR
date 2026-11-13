# Operator overloading for opt_expr.

new_relation <- function(lhs, rhs, sense) {
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
#' closed feasible sets. Logical operators are not overloaded; R's `&&` and
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
#' @param ... Unused.
#' @export
Math.opt_expr <- function(x, ...) {
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
#' @param y An expression or numeric matrix.
#' @export
matrixOps.opt_expr <- function(x, y) {
  if (.Generic != "%*%") stop("Unsupported matrix operator.", call. = FALSE)
  expr_matmul(x, y)
}

as_matrix_shape <- function(e, left) {
  if (!is.null(e$dim) && length(e$dim) == 2L) return(e$dim)
  if (left) c(1L, e$len) else c(e$len, 1L)
}

expr_matmul <- function(x, y) {
  if (!is_expr(x) && !is_expr(y)) return(x %*% y)
  if (is_expr(x) && is_expr(y)) {
    dx <- as_matrix_shape(x, TRUE)
    dy <- as_matrix_shape(y, FALSE)
    if (dx[2L] != dy[1L]) stop("Non-conformable arguments in %*%.",
                               call. = FALSE)
    if (dx[1L] == 1L && dy[2L] == 1L) return(expr_sum(expr_mul(x, y)))
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
    A <- x
    if (is.null(dim(A))) A <- matrix(A, nrow = 1L)
    dy <- as_matrix_shape(y, FALSE)
    if (ncol(A) != dy[1L]) stop("Non-conformable arguments in %*%.",
                                call. = FALSE)
    M <- Matrix::kronecker(Matrix::Diagonal(dy[2L]), Matrix::Matrix(A, sparse = TRUE, doDiag = FALSE))
    return(expr_linmap(M, y, dim = c(nrow(A), dy[2L])))
  }
  B <- y
  if (is.null(dim(B))) B <- matrix(B, ncol = 1L)
  dx <- as_matrix_shape(x, TRUE)
  if (dx[2L] != nrow(B)) stop("Non-conformable arguments in %*%.",
                              call. = FALSE)
  M <- Matrix::kronecker(Matrix::t(Matrix::Matrix(B, sparse = TRUE, doDiag = FALSE)),
                         Matrix::Diagonal(dx[1L]))
  expr_linmap(M, x, dim = c(dx[1L], ncol(B)))
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
  if (nidx == 1L && is.numeric(i) && length(i) && !anyNA(i) &&
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
  out <- NULL
  for (k in seq_along(args)) {
    a <- args[[k]]
    if (a$type == "poly") {
      a <- new_poly(a$model, n, NULL,
                    const = replace(numeric(n), off[k] + seq_len(a$len), a$const),
                    li = a$li + off[k], lj = a$lj, lv = a$lv,
                    qi = a$qi + off[k], qj = a$qj, qk = a$qk, qv = a$qv)
    } else {
      M <- Matrix::sparseMatrix(i = off[k] + seq_len(a$len),
                                j = seq_len(a$len), x = 1,
                                dims = c(n, a$len))
      a <- expr_linmap(M, a)
    }
    out <- if (is.null(out)) a else expr_add(out, a)
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
