# Expression representation.
#
# Every opt_expr is a vector of `len` scalar elements (with an optional dim).
# Two node types exist:
#
# * "poly": a polynomial of degree at most two in the model's symbols
#   (decision variables and parameters share one symbol space). Stored as
#   triplets: constant vector `const`, linear terms (li, lj, lv) meaning
#   element li gets lv * z[lj], and quadratic terms (qi, qj, qk, qv) meaning
#   element qi gets qv * z[qj] * z[qk]. Duplicate triplets add; compaction is
#   lazy.
# * "nl": a nonlinear node with an operator `op`, a list of argument
#   expressions `args`, and operator data `data`.
#
# Expressions are immutable values. They hold a reference to their model
# environment (or NULL for pure constants).

new_poly <- function(model, len, dim = NULL, const = numeric(len),
                     li = integer(), lj = integer(), lv = numeric(),
                     qi = integer(), qj = integer(), qk = integer(),
                     qv = numeric()) {
  # class<- instead of structure(): this constructor is on the hot path
  e <- list(type = "poly", model = model, len = as.integer(len),
            dim = dim, const = const, li = li, lj = lj, lv = lv,
            qi = qi, qj = qj, qk = qk, qv = qv)
  class(e) <- "opt_expr"
  e
}

new_nl <- function(model, len, dim, op, args, data = NULL) {
  structure(list(type = "nl", model = model, len = as.integer(len),
                 dim = dim, op = op, args = args, data = data),
            class = "opt_expr")
}

is_expr <- function(x) inherits(x, "opt_expr")

#' Test whether an object is an optimization expression
#' @param x Any object.
#' @return `TRUE` if `x` is an `opt_expr`.
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x")
#' is_opt_expr(2 * x + 1)
#' is_opt_expr(3)
#' @export
is_opt_expr <- function(x) is_expr(x)

# Coerce numeric input to a constant poly expression.
as_expr <- function(x) {
  if (is_expr(x)) return(x)
  if (inherits(x, "Matrix")) x <- as.matrix(x)
  if (is.logical(x)) x <- as.numeric(x)
  if (!is.numeric(x)) {
    stop("Cannot combine an optimization expression with an object of class ",
         paste(class(x), collapse = "/"), ".", call. = FALSE)
  }
  if (anyNA(x)) stop("Missing values are not allowed in expressions.",
                     call. = FALSE)
  d <- dim(x)
  if (!is.null(d) && length(d) == 1L) d <- NULL
  new_poly(NULL, length(x), d, const = as.vector(x))
}

join_models <- function(a, b) {
  ma <- a$model
  mb <- b$model
  if (is.null(ma)) return(mb)
  if (is.null(mb) || identical(ma, mb)) return(ma)
  stop("Cannot combine expressions that belong to different models.",
       call. = FALSE)
}

join_len_dim <- function(a, b) {
  la <- a$len
  lb <- b$len
  if (la != lb && la != 1L && lb != 1L) {
    stop(sprintf(
      "Non-conformable expressions: lengths %d and %d (only length one is recycled).",
      la, lb), call. = FALSE)
  }
  n <- max(la, lb)
  da <- if (la == n) a$dim else NULL
  db <- if (lb == n) b$dim else NULL
  if (!is.null(da) && !is.null(db) && !identical(as.integer(da), as.integer(db))) {
    stop("Non-conformable expressions: dimensions ",
         paste(da, collapse = "x"), " and ", paste(db, collapse = "x"), ".",
         call. = FALSE)
  }
  list(len = n, dim = da %||% db)
}

poly_degree <- function(e) {
  if (length(e$qv)) 2L else if (length(e$lv)) 1L else 0L
}

is_constant_poly <- function(e) {
  e$type == "poly" && length(e$lv) == 0L && length(e$qv) == 0L
}

# Merge duplicate triplets and drop exact zeros.
compact_poly <- function(e) {
  if (e$type != "poly") return(e)
  if (length(e$lv)) {
    key <- paste(e$li, e$lj)
    if (anyDuplicated(key) || any(e$lv == 0)) {
      s <- rowsum(e$lv, key, reorder = FALSE)
      first <- !duplicated(key)
      keep <- s[, 1L] != 0
      e$li <- e$li[first][keep]
      e$lj <- e$lj[first][keep]
      e$lv <- unname(s[keep, 1L])
    }
  }
  if (length(e$qv)) {
    j <- pmin(e$qj, e$qk)
    k <- pmax(e$qj, e$qk)
    key <- paste(e$qi, j, k)
    s <- rowsum(e$qv, key, reorder = FALSE)
    first <- !duplicated(key)
    keep <- s[, 1L] != 0
    e$qi <- e$qi[first][keep]
    e$qj <- j[first][keep]
    e$qk <- k[first][keep]
    e$qv <- unname(s[keep, 1L])
  }
  e
}

# Select (and possibly repeat) elements of an expression.
expr_subset <- function(e, idx, dim = NULL) {
  idx <- as.integer(idx)
  if (anyNA(idx) || any(idx < 1L | idx > e$len)) {
    stop("Expression index out of bounds.", call. = FALSE)
  }
  n <- length(idx)
  if (e$type == "nl") {
    if (n == e$len && identical(idx, seq_len(e$len))) {
      e$dim <- dim
      return(e)
    }
    M <- Matrix::sparseMatrix(i = seq_len(n), j = idx, x = 1,
                              dims = c(n, e$len))
    return(new_nl(e$model, n, dim, "linmap", list(e), list(M = M)))
  }
  gl <- gather_index(e$li, e$len, idx)
  gq <- gather_index(e$qi, e$len, idx)
  new_poly(e$model, n, dim, const = e$const[idx],
           li = gl$i, lj = e$lj[gl$pos], lv = e$lv[gl$pos],
           qi = gq$i, qj = e$qj[gq$pos], qk = e$qk[gq$pos], qv = e$qv[gq$pos])
}

recycle_expr <- function(e, n) {
  if (e$len == n) return(e)
  expr_subset(e, rep.int(1L, n))
}

# Multiply each element by a constant (vector of length 1 or len).
expr_scale <- function(e, s) {
  if (length(s) != 1L && length(s) != e$len) {
    stop("Internal error: scale length mismatch.", call. = FALSE)
  }
  if (e$type == "nl") {
    cst <- new_poly(NULL, length(s), NULL, const = s)
    return(nl_binary("*", cst, e))
  }
  if (length(s) == 1L) {
    e$const <- e$const * s
    e$lv <- e$lv * s
    e$qv <- e$qv * s
  } else {
    e$const <- e$const * s
    e$lv <- e$lv * s[e$li]
    e$qv <- e$qv * s[e$qi]
  }
  e
}

expr_add <- function(a, b, sign = 1) {
  a <- as_expr(a)
  b <- as_expr(b)
  model <- join_models(a, b)
  ld <- join_len_dim(a, b)
  if (a$type == "nl" || b$type == "nl") {
    if (sign < 0) b <- expr_scale(b, -1)
    return(nl_binary("+", a, b))
  }
  n <- ld$len
  a <- recycle_expr(a, n)
  b <- recycle_expr(b, n)
  new_poly(model, n, ld$dim,
           const = a$const + sign * b$const,
           li = c(a$li, b$li), lj = c(a$lj, b$lj), lv = c(a$lv, sign * b$lv),
           qi = c(a$qi, b$qi), qj = c(a$qj, b$qj), qk = c(a$qk, b$qk),
           qv = c(a$qv, sign * b$qv))
}

# Product of two poly expressions of total degree <= 2 (same length).
poly_product <- function(a, b, model, n, dim) {
  # constant parts
  const <- a$const * b$const
  li <- c(a$li, b$li)
  lj <- c(a$lj, b$lj)
  lv <- c(a$lv * b$const[a$li], b$lv * a$const[b$li])
  qi <- c(a$qi, b$qi)
  qj <- c(a$qj, b$qj)
  qk <- c(a$qk, b$qk)
  qv <- c(a$qv * b$const[a$qi], b$qv * a$const[b$qi])
  if (length(a$lv) && length(b$lv)) {
    ca <- tabulate(a$li, nbins = n)
    cb <- tabulate(b$li, nbins = n)
    oa <- order(a$li)
    ob <- order(b$li)
    sa <- cumsum(c(0L, ca))[seq_len(n)]
    sb <- cumsum(c(0L, cb))[seq_len(n)]
    np <- ca * cb
    if (sum(np) > 0L) {
      el <- rep.int(seq_len(n), np)
      t <- sequence(np) - 1L
      ai <- oa[sa[el] + t %/% cb[el] + 1L]
      bi <- ob[sb[el] + t %% cb[el] + 1L]
      qi <- c(qi, el)
      qj <- c(qj, a$lj[ai])
      qk <- c(qk, b$lj[bi])
      qv <- c(qv, a$lv[ai] * b$lv[bi])
    }
  }
  new_poly(model, n, dim, const = const, li = li, lj = lj, lv = lv,
           qi = qi, qj = qj, qk = qk, qv = qv)
}

expr_mul <- function(a, b) {
  a <- as_expr(a)
  b <- as_expr(b)
  model <- join_models(a, b)
  ld <- join_len_dim(a, b)
  n <- ld$len
  if (is_constant_poly(a) || is_constant_poly(b)) {
    if (is_constant_poly(b)) {
      tmp <- a
      a <- b
      b <- tmp
    }
    b <- recycle_expr(b, n)
    out <- expr_scale(b, if (a$len == 1L) a$const else rep_len(a$const, n))
    out$dim <- ld$dim
    out$model <- model
    return(out)
  }
  if (a$type == "poly" && b$type == "poly") {
    a <- compact_poly(a)
    b <- compact_poly(b)
    if (poly_degree(a) + poly_degree(b) <= 2L) {
      return(poly_product(recycle_expr(a, n), recycle_expr(b, n), model, n,
                          ld$dim))
    }
  }
  nl_binary("*", a, b)
}

expr_div <- function(a, b) {
  a <- as_expr(a)
  b <- as_expr(b)
  if (is_constant_poly(b)) {
    if (any(b$const == 0)) {
      stop("Division of an expression by zero.", call. = FALSE)
    }
    return(expr_mul(a, new_poly(NULL, b$len, b$dim, const = 1 / b$const)))
  }
  nl_binary("/", a, b)
}

expr_pow <- function(a, p) {
  a <- as_expr(a)
  if (is_expr(p) && !is_constant_poly(p)) {
    # a^p with a variable exponent is exp(p * log(a)), defined for a > 0
    return(expr_unary("exp", expr_mul(p, expr_unary("log", a))))
  }
  pv <- if (is_expr(p)) p$const else p
  if (!is.numeric(pv) || anyNA(pv)) {
    stop("Exponent must be numeric.", call. = FALSE)
  }
  if (is_constant_poly(a)) {
    pe <- as_expr(pv)
    ld <- join_len_dim(a, pe)
    return(new_poly(NULL, ld$len, ld$dim,
                    const = rep_len(a$const, ld$len)^rep_len(pv, ld$len)))
  }
  if (length(pv) == 1L) {
    if (pv == 0) {
      return(new_poly(a$model, a$len, a$dim, const = rep(1, a$len)))
    }
    if (pv == 1) return(a)
    if (pv == 2 && a$type == "poly") {
      a <- compact_poly(a)
      if (poly_degree(a) <= 1L) return(poly_product(a, a, a$model, a$len,
                                                    a$dim))
    }
  }
  expr_pow_const(a, pv)
}

expr_sum <- function(e) {
  if (e$type == "nl") {
    M <- Matrix::sparseMatrix(i = rep.int(1L, e$len), j = seq_len(e$len),
                              x = 1, dims = c(1L, e$len))
    return(new_nl(e$model, 1L, NULL, "linmap", list(e), list(M = M)))
  }
  new_poly(e$model, 1L, NULL, const = sum(e$const),
           li = rep.int(1L, length(e$li)), lj = e$lj, lv = e$lv,
           qi = rep.int(1L, length(e$qi)), qj = e$qj, qk = e$qk, qv = e$qv)
}

# Apply a constant linear map M (k x len) to an expression.
expr_linmap <- function(M, e, dim = NULL) {
  M <- methods::as(Matrix::Matrix(M, sparse = TRUE), "CsparseMatrix")
  if (ncol(M) != e$len) {
    stop(sprintf("Non-conformable: matrix has %d columns, expression has %d elements.",
                 ncol(M), e$len), call. = FALSE)
  }
  k <- nrow(M)
  if (e$type == "nl") return(new_nl(e$model, k, dim, "linmap", list(e),
                                    list(M = M)))
  const <- as.vector(M %*% e$const)
  out <- new_poly(e$model, k, dim, const = const)
  ns <- max(c(e$lj, e$qj, e$qk, 1L))
  if (length(e$lv)) {
    L <- Matrix::sparseMatrix(i = e$li, j = e$lj, x = e$lv,
                              dims = c(e$len, ns))
    T <- methods::as(M %*% L, "TsparseMatrix")
    out$li <- T@i + 1L
    out$lj <- T@j + 1L
    out$lv <- T@x
  }
  if (length(e$qv)) {
    j <- pmin(e$qj, e$qk)
    kk <- pmax(e$qj, e$qk)
    key <- (kk - 1) * ns + j
    ukey <- unique(key)
    Q <- Matrix::sparseMatrix(i = e$qi, j = match(key, ukey), x = e$qv,
                              dims = c(e$len, length(ukey)))
    T <- methods::as(M %*% Q, "TsparseMatrix")
    kv <- ukey[T@j + 1L]
    out$qi <- T@i + 1L
    out$qj <- as.integer((kv - 1) %% ns + 1)
    out$qk <- as.integer((kv - 1) %/% ns + 1)
    out$qv <- T@x
  }
  out
}

# Polynomial expressions that are linear in the symbols.
expr_kind_label <- function(e) {
  if (e$type == "nl") return("nonlinear")
  switch(as.character(poly_degree(compact_poly(e))),
         "0" = "constant", "1" = "affine", "2" = "quadratic")
}
