# Nonlinear expression nodes: construction, evaluation, exact first and
# second derivatives (forward mode sparse Jacobians and weighted Hessians),
# and parameter lowering.

unary_ops <- c("exp", "log", "sqrt", "sin", "cos", "tan", "tanh", "sinh",
               "cosh", "abs", "log1p", "expm1", "atan", "asin", "acos")

nl_binary <- function(op, a, b) {
  model <- join_models(a, b)
  ld <- join_len_dim(a, b)
  new_nl(model, ld$len, ld$dim, op, list(a, b))
}

expr_unary <- function(op, a) {
  a <- as_expr(a)
  if (is_constant_poly(a)) {
    return(new_poly(NULL, a$len, a$dim, const = match.fun(op)(a$const)))
  }
  new_nl(a$model, a$len, a$dim, op, list(a))
}

# Constant exponent power node (general exponents are rewritten through
# exp and log by expr_pow).
expr_pow_const <- function(a, p) {
  ld <- join_len_dim(a, as_expr(p))
  new_nl(a$model, ld$len, ld$dim, "^", list(a), list(p = p))
}

unary_d1 <- function(op, g) {
  switch(op,
    exp = exp(g), log = 1 / g, sqrt = 0.5 / sqrt(g), sin = cos(g),
    cos = -sin(g), tan = 1 / cos(g)^2, tanh = 1 - tanh(g)^2,
    sinh = cosh(g), cosh = sinh(g), abs = sign(g), log1p = 1 / (1 + g),
    expm1 = exp(g), atan = 1 / (1 + g^2), asin = 1 / sqrt(1 - g^2),
    acos = -1 / sqrt(1 - g^2))
}

unary_d2 <- function(op, g) {
  switch(op,
    exp = exp(g), log = -1 / g^2, sqrt = -0.25 / g^1.5,
    sin = -sin(g), cos = -cos(g), tan = 2 * tan(g) / cos(g)^2,
    tanh = -2 * tanh(g) * (1 - tanh(g)^2), sinh = sinh(g), cosh = cosh(g),
    abs = 0 * g, log1p = -1 / (1 + g)^2, expm1 = exp(g),
    atan = -2 * g / (1 + g^2)^2, asin = g / (1 - g^2)^1.5,
    acos = -g / (1 - g^2)^1.5)
}

# ---------------------------------------------------------------- evaluation

eval_expr <- function(e, z) {
  if (e$type == "poly") {
    out <- e$const
    if (length(e$lv)) out <- out + acc_sum(e$li, e$lv * z[e$lj], e$len)
    if (length(e$qv)) {
      out <- out + acc_sum(e$qi, e$qv * z[e$qj] * z[e$qk], e$len)
    }
    return(out)
  }
  args <- lapply(e$args, eval_expr, z = z)
  n <- e$len
  switch(e$op,
    "+" = rep_len(args[[1L]] + args[[2L]], n),
    "*" = rep_len(args[[1L]] * args[[2L]], n),
    "/" = rep_len(args[[1L]] / args[[2L]], n),
    "^" = rep_len(args[[1L]]^e$data$p, n),
    linmap = as.vector(e$data$M %*% args[[1L]]),
    fn = e$data$f(args[[1L]]),
    match.fun(e$op)(args[[1L]]))
}

rep_rows <- function(J, n) {
  if (nrow(J) == n) return(J)
  J[rep.int(1L, n), , drop = FALSE]
}

diag_mul <- function(d, J) {
  Matrix::Diagonal(x = d) %*% J
}

poly_jac <- function(e, z, ns) {
  i <- c(e$li, e$qi, e$qi)
  j <- c(e$lj, e$qj, e$qk)
  x <- c(e$lv, e$qv * z[e$qk], e$qv * z[e$qj])
  Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(e$len, ns))
}

# Value and sparse Jacobian (len x ns) of an expression at z.
jac_expr <- function(e, z, ns) {
  if (e$type == "poly") {
    return(list(val = eval_expr(e, z), J = poly_jac(e, z, ns)))
  }
  n <- e$len
  if (e$op == "linmap") {
    a <- jac_expr(e$args[[1L]], z, ns)
    M <- e$data$M
    return(list(val = as.vector(M %*% a$val), J = M %*% a$J))
  }
  if (e$op == "fn") {
    a <- jac_expr(e$args[[1L]], z, ns)
    g <- fn_node_grad(e, a$val)
    return(list(val = e$data$f(a$val), J = Matrix::Matrix(g, 1L) %*% a$J))
  }
  if (e$op %in% unary_ops) {
    a <- jac_expr(e$args[[1L]], z, ns)
    return(list(val = match.fun(e$op)(a$val),
                J = diag_mul(unary_d1(e$op, a$val), a$J)))
  }
  if (e$op == "^") {
    a <- jac_expr(e$args[[1L]], z, ns)
    p <- e$data$p
    g <- rep_len(a$val, n)
    p <- rep_len(p, n)
    return(list(val = g^p, J = diag_mul(pow_d1(g, p), rep_rows(a$J, n))))
  }
  a <- jac_expr(e$args[[1L]], z, ns)
  b <- jac_expr(e$args[[2L]], z, ns)
  av <- rep_len(a$val, n)
  bv <- rep_len(b$val, n)
  Ja <- rep_rows(a$J, n)
  Jb <- rep_rows(b$J, n)
  switch(e$op,
    "+" = list(val = av + bv, J = Ja + Jb),
    "*" = list(val = av * bv, J = diag_mul(bv, Ja) + diag_mul(av, Jb)),
    "/" = list(val = av / bv,
               J = diag_mul(1 / bv, Ja) - diag_mul(av / bv^2, Jb)))
}

# Weighted Hessian sum_e w[e] * Hessian(e-th element) at z (ns x ns).
hess_expr <- function(e, z, ns, w) {
  w <- rep_len(w, e$len)
  if (e$type == "poly") {
    if (!length(e$qv)) return(Matrix::sparseMatrix(i = integer(), j = integer(),
                                                   x = numeric(),
                                                   dims = c(ns, ns)))
    v <- e$qv * w[e$qi]
    return(Matrix::sparseMatrix(i = c(e$qj, e$qk), j = c(e$qk, e$qj),
                                x = c(v, v), dims = c(ns, ns)))
  }
  n <- e$len
  agg <- function(arg, wv) {
    # weights for an argument that may be recycled from length one
    if (arg$len == 1L && n > 1L) sum(wv) else wv
  }
  cross <- function(Ja, Jb, d) {
    Ja <- rep_rows(Ja, n)
    Jb <- rep_rows(Jb, n)
    Matrix::t(Ja) %*% diag_mul(d, Jb)
  }
  if (e$op == "linmap") {
    return(hess_expr(e$args[[1L]], z, ns,
                     as.vector(Matrix::crossprod(e$data$M, w))))
  }
  if (e$op == "fn") {
    if (is.null(e$data$hess)) {
      stop(sprintf("No Hessian was supplied for external function '%s'.",
                   e$data$name), call. = FALSE)
    }
    arg <- e$args[[1L]]
    a <- jac_expr(arg, z, ns)
    Hf <- Matrix::Matrix(as.matrix(e$data$hess(a$val)) * w)
    return(Matrix::t(a$J) %*% Hf %*% a$J +
             hess_expr(arg, z, ns, w * fn_node_grad(e, a$val)))
  }
  if (e$op %in% unary_ops || e$op == "^") {
    arg <- e$args[[1L]]
    a <- jac_expr(arg, z, ns)
    g <- rep_len(a$val, n)
    if (e$op == "^") {
      p <- rep_len(e$data$p, n)
      d1 <- pow_d1(g, p)
      d2 <- pow_d2(g, p)
    } else {
      d1 <- unary_d1(e$op, g)
      d2 <- unary_d2(e$op, g)
    }
    return(cross(a$J, a$J, w * d2) + hess_expr(arg, z, ns, agg(arg, w * d1)))
  }
  A <- e$args[[1L]]
  B <- e$args[[2L]]
  if (e$op == "+") {
    return(hess_expr(A, z, ns, agg(A, w)) + hess_expr(B, z, ns, agg(B, w)))
  }
  a <- jac_expr(A, z, ns)
  b <- jac_expr(B, z, ns)
  av <- rep_len(a$val, n)
  bv <- rep_len(b$val, n)
  if (e$op == "*") {
    H <- hess_expr(A, z, ns, agg(A, w * bv)) +
      hess_expr(B, z, ns, agg(B, w * av)) +
      cross(a$J, b$J, w) + cross(b$J, a$J, w)
    return(H)
  }
  # division a / b
  H <- hess_expr(A, z, ns, agg(A, w / bv)) +
    hess_expr(B, z, ns, agg(B, -w * av / bv^2)) -
    cross(a$J, b$J, w / bv^2) - cross(b$J, a$J, w / bv^2) +
    cross(b$J, b$J, 2 * w * av / bv^3)
  H
}

# ------------------------------------------------------------------ lowering

# Substitute current parameter values into an expression. The result only
# depends on decision variables; nonlinear nodes are rebuilt through the
# ordinary constructors, so a node that becomes polynomial after
# substitution (for example p * x^2 with p fixed) collapses to a poly node.
lower_expr <- function(e, is_par, pval) {
  if (e$type == "poly") {
    if (!any(is_par[c(e$lj, e$qj, e$qk)])) return(e)
    const <- e$const
    lp <- is_par[e$lj]
    if (any(lp)) {
      const <- const + acc_sum(e$li[lp], e$lv[lp] * pval[e$lj[lp]], e$len)
    }
    pj <- is_par[e$qj]
    pk <- is_par[e$qk]
    both <- pj & pk
    if (any(both)) {
      const <- const + acc_sum(e$qi[both],
                               e$qv[both] * pval[e$qj[both]] * pval[e$qk[both]],
                               e$len)
    }
    onej <- pj & !pk
    onek <- pk & !pj
    li <- c(e$li[!lp], e$qi[onej], e$qi[onek])
    lj <- c(e$lj[!lp], e$qk[onej], e$qj[onek])
    lv <- c(e$lv[!lp], e$qv[onej] * pval[e$qj[onej]],
            e$qv[onek] * pval[e$qk[onek]])
    keep <- !(pj | pk)
    return(new_poly(e$model, e$len, e$dim, const = const, li = li, lj = lj,
                    lv = lv, qi = e$qi[keep], qj = e$qj[keep],
                    qk = e$qk[keep], qv = e$qv[keep]))
  }
  args <- lapply(e$args, lower_expr, is_par = is_par, pval = pval)
  out <- switch(e$op,
    "+" = expr_add(args[[1L]], args[[2L]]),
    "*" = expr_mul(args[[1L]], args[[2L]]),
    "/" = expr_div(args[[1L]], args[[2L]]),
    "^" = expr_pow(args[[1L]], e$data$p),
    linmap = expr_linmap(e$data$M, args[[1L]]),
    fn = new_nl(e$model, 1L, NULL, "fn", args, e$data),
    expr_unary(e$op, args[[1L]]))
  out$dim <- e$dim
  out
}

# Symbols used by an expression.
expr_symbols <- function(e) {
  if (e$type == "poly") return(sort(unique(c(e$lj, e$qj, e$qk))))
  sort(unique(unlist(lapply(e$args, expr_symbols), use.names = FALSE)))
}

expr_has_nl <- function(e) e$type == "nl"

# Derivatives of g^p, elementwise in p, exact for p in {0, 1, 2} even at
# g = 0 where the general formula would evaluate 0 * 0^-1.
pow_d1 <- function(g, p) {
  ifelse(p == 0, 0, ifelse(p == 1, 1, p * g^(p - 1)))
}
pow_d2 <- function(g, p) {
  ifelse(p == 0 | p == 1, 0, ifelse(p == 2, 2, p * (p - 1) * g^(p - 2)))
}
