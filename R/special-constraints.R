# Constraints that are not rows l <= g(x) <= u: cones, indicator
# constraints, special ordered sets and complementarity. They are stored in
# model$special and lowered into the standard form next to the rows.

add_special <- function(model, kind, item, name) {
  sp <- model$special %||% list()
  have <- vapply(sp[[kind]] %||% list(), `[[`, "", "name")
  if (is.null(name)) {
    k <- length(have) + 1L
    repeat {
      name <- paste0(kind, k)
      if (!name %in% have) break
      k <- k + 1L
    }
  }
  check_scalar_string(name, "name")
  if (name %in% have) {
    stop(sprintf("The model already has a %s constraint named '%s'.", kind,
                 name), call. = FALSE)
  }
  item$name <- name
  sp[[kind]] <- c(sp[[kind]], list(item))
  model$special <- sp
  touch_structure(model)
  invisible(structure(list(model = model, kind = kind, name = name),
                      class = "opt_special_constraint"))
}

#' @export
print.opt_special_constraint <- function(x, ...) {
  cat(sprintf("<opt_special_constraint %s '%s'>\n", x$kind, x$name))
  invisible(x)
}

model_expr <- function(model, x, what) {
  x <- as_expr(x)
  if (!is.null(x$model) && !identical(x$model, model)) {
    stop(sprintf("`%s` belongs to a different model.", what), call. = FALSE)
  }
  x
}

check_affine <- function(x, what) {
  if (x$type != "poly" || length(x$qv)) {
    stop(sprintf("`%s` must be affine (linear plus a constant).", what),
         call. = FALSE)
  }
  x
}

# ------------------------------------------------------------------ cones

#' Conic constraints
#'
#' Convex cone constraints on affine expressions, solved by conic solvers
#' (Clarabel) and, for second-order cones, by SCIP:
#'
#' * `opt_soc_constraint()`: \eqn{\|x\|_2 \le t} (second-order cone).
#' * `opt_exp_cone()`: \eqn{y \exp(x / y) \le z} with \eqn{y > 0}
#'   (the closure of the exponential cone, as in Clarabel and SCS).
#' * `opt_pow_cone()`: \eqn{x^\alpha y^{1-\alpha} \ge |z|} with
#'   \eqn{x, y \ge 0} and \eqn{0 < \alpha < 1}.
#' * `opt_psd_constraint()`: the symmetric matrix `M` is positive
#'   semidefinite. Only the upper triangle of `M` is used, so symmetry must
#'   hold by construction (for example with [opt_sym_variable()]).
#'
#' All arguments must be affine in the decision variables after parameter
#' substitution. Models with cones are classified as `"SOCP"` (second-order
#' cones only) or `"CONIC"`, with `"MI"` prefixed when integer variables are
#' present; they are convex when the objective is linear or a convex
#' quadratic and all rows are linear.
#'
#' @param model An [opt_model()].
#' @param x,y,z,t Affine expressions; `x` in `opt_soc_constraint()` may be a
#'   vector, all others are scalars.
#' @param alpha Exponent of the power cone.
#' @param M Square affine matrix expression.
#' @param name Optional name.
#' @return An `opt_special_constraint` reference, invisibly.
#' @examplesIf requireNamespace("clarabel", quietly = TRUE)
#' # distance from the origin to the line x + y = 2
#' m <- opt_model()
#' x <- opt_variable(m, "x", dim = 2)
#' t <- opt_variable(m, "t")
#' opt_constraint(m, sum(x) == 2)
#' opt_soc_constraint(m, x, t)
#' opt_objective(m, t)
#' opt_objective_value(opt_solve(m))  # sqrt(2)
#' @export
opt_soc_constraint <- function(model, x, t, name = NULL) {
  check_model(model)
  x <- check_affine(model_expr(model, x, "x"), "x")
  t <- check_affine(model_expr(model, t, "t"), "t")
  if (t$len != 1L) stop("`t` must be scalar.", call. = FALSE)
  add_special(model, "cone", list(type = "soc", args = list(t, x)), name)
}

#' @rdname opt_soc_constraint
#' @export
opt_exp_cone <- function(model, x, y, z, name = NULL) {
  check_model(model)
  args <- lapply(list(x = x, y = y, z = z), function(a) a)
  args <- Map(function(a, nm) {
    a <- check_affine(model_expr(model, a, nm), nm)
    if (a$len != 1L) stop(sprintf("`%s` must be scalar.", nm), call. = FALSE)
    a
  }, args, names(args))
  add_special(model, "cone", list(type = "exp", args = unname(args)), name)
}

#' @rdname opt_soc_constraint
#' @export
opt_pow_cone <- function(model, x, y, z, alpha, name = NULL) {
  check_model(model)
  check_numeric(alpha, "alpha")
  if (length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a number strictly between 0 and 1.", call. = FALSE)
  }
  args <- Map(function(a, nm) {
    a <- check_affine(model_expr(model, a, nm), nm)
    if (a$len != 1L) stop(sprintf("`%s` must be scalar.", nm), call. = FALSE)
    a
  }, list(x, y, z), c("x", "y", "z"))
  add_special(model, "cone", list(type = "pow", args = args, alpha = alpha),
              name)
}

#' @rdname opt_soc_constraint
#' @export
opt_psd_constraint <- function(model, M, name = NULL) {
  check_model(model)
  M <- check_affine(model_expr(model, M, "M"), "M")
  d <- M$dim
  if (length(d) != 2L || d[1L] != d[2L]) {
    stop("`M` must be a square matrix expression.", call. = FALSE)
  }
  k <- d[1L]
  # upper triangle by columns, as used by the conic solvers
  idx <- which(upper.tri(matrix(0, k, k), diag = TRUE))
  vec <- expr_subset(M, idx, NULL)
  add_special(model, "cone", list(type = "psd", args = list(vec), k = k),
              name)
}

#' Symmetric matrix variables
#'
#' Creates \eqn{k(k+1)/2} decision variables and returns them as a
#' symmetric \eqn{k \times k} matrix expression, for use with
#' [opt_psd_constraint()] and matrix products.
#'
#' @inheritParams opt_variable
#' @param k Matrix size.
#' @return A `k` by `k` `opt_expr`. The underlying variable block is named
#'   `name` and holds the upper triangle by columns.
#' @examplesIf requireNamespace("clarabel", quietly = TRUE)
#' # smallest t with t * I - A positive semidefinite: the largest eigenvalue
#' A <- matrix(c(2, 1, 1, 3), 2)
#' m <- opt_model()
#' X <- opt_sym_variable(m, "X", 2)
#' t <- opt_variable(m, "t")
#' opt_constraint(m, X == t * diag(2) - A)
#' opt_psd_constraint(m, X)
#' opt_objective(m, t)
#' opt_objective_value(opt_solve(m))
#' max(eigen(A)$values)
#' @export
opt_sym_variable <- function(model, name, k, lower = NULL, upper = NULL) {
  check_model(model)
  if (!is.numeric(k) || length(k) != 1L || k < 1 || k != round(k)) {
    stop("`k` must be a positive integer.", call. = FALSE)
  }
  k <- as.integer(k)
  v <- opt_variable(model, name, dim = k * (k + 1L) / 2L, lower = lower,
                    upper = upper)
  pos <- matrix(0L, k, k)
  pos[upper.tri(pos, diag = TRUE)] <- seq_len(k * (k + 1L) / 2L)
  pos[lower.tri(pos)] <- t(pos)[lower.tri(pos)]
  out <- expr_subset(v, as.vector(pos), NULL)
  out$dim <- c(k, k)
  out
}

# -------------------------------------------------------- discrete forms

as_binary_symbol <- function(model, z, what) {
  sym <- resolve_symbols(model, z, "var")
  if (!all(model$type[sym] == "binary")) {
    stop(sprintf("`%s` must be binary variables.", what), call. = FALSE)
  }
  sym
}

#' Indicator constraints
#'
#' `opt_indicator(model, z, relation)` imposes the linear `relation` only
#' when the binary variable `z` equals `active`; otherwise the relation is
#' free. SCIP handles indicator constraints natively. For other
#' mixed-integer solvers they are reformulated exactly with big-M rows,
#' where M is computed from the variable bounds; that requires finite bounds
#' on the variables of the relation, and the solve stops with an error
#' otherwise. No approximation is made.
#'
#' @param model An [opt_model()].
#' @param z A binary variable (scalar) or a vector of binaries, one per
#'   element of `relation`.
#' @param relation A linear relation such as `x + y <= 4`.
#' @param active `1` (the relation holds when `z = 1`) or `0`.
#' @param name Optional name.
#' @return An `opt_special_constraint` reference, invisibly.
#' @examples
#' # a fixed cost of 10 is paid when production x is positive
#' m <- opt_model()
#' x <- opt_variable(m, "x", lower = 0, upper = 8)
#' on <- opt_variable(m, "on", type = "binary")
#' opt_indicator(m, on, x <= 0, active = 0)
#' opt_constraint(m, x >= 3)
#' opt_objective(m, 2 * x + 10 * on)
#' opt_value(opt_solve(m), c(x, on))
#' @export
opt_indicator <- function(model, z, relation, active = 1, name = NULL) {
  check_model(model)
  if (!inherits(relation, "opt_relation")) {
    stop("`relation` must be a comparison such as `x <= 4`.", call. = FALSE)
  }
  zs <- as_binary_symbol(model, z, "z")
  e <- check_affine(model_expr(model, relation$expr, "relation"), "relation")
  if (!active %in% c(0, 1)) stop("`active` must be 0 or 1.", call. = FALSE)
  if (length(zs) != 1L && length(zs) != e$len) {
    stop("`z` must be scalar or have one element per relation element.",
         call. = FALSE)
  }
  add_special(model, "indicator",
              list(z = rep_len(zs, e$len), expr = e, sense = relation$sense,
                   active = active), name)
}

#' Special ordered sets and complementarity
#'
#' `opt_sos1()`: at most one element of `x` is nonzero. `opt_sos2()`: at
#' most two elements are nonzero and they are adjacent in the order given by
#' `weights`. `opt_complementarity()`: \eqn{x_i y_i = 0} elementwise for
#' nonnegative variables, the complementarity condition
#' \eqn{0 \le x \perp y \ge 0}, as SOS1 pairs.
#'
#' SCIP handles these natively. For other mixed-integer solvers they are
#' reformulated exactly with one binary per element (SOS1) or segment
#' (SOS2), which requires finite bounds.
#'
#' @param model An [opt_model()].
#' @param x,y Decision variables (a block or a slice of one).
#' @param weights Distinct weights that order the set; defaults to
#'   `seq_along(x)`.
#' @param name Optional name.
#' @return An `opt_special_constraint` reference (a list of them for
#'   `opt_complementarity()`), invisibly.
#' @examples
#' # piecewise linear interpolation of f(u) = u^2 on the grid 0, 1, 2, 3
#' m <- opt_model()
#' lam <- opt_variable(m, "lambda", dim = 4, lower = 0, upper = 1)
#' u <- opt_variable(m, "u")
#' opt_constraint(m, sum(lam) == 1)
#' opt_constraint(m, u == sum(0:3 * lam))
#' opt_constraint(m, u >= 1.5)
#' opt_sos2(m, lam)
#' opt_objective(m, sum((0:3)^2 * lam))
#' opt_objective_value(opt_solve(m))  # 2.5, the interpolant at 1.5
#' @export
opt_sos1 <- function(model, x, weights = NULL, name = NULL) {
  sos_add(model, x, weights, name, 1L)
}

#' @rdname opt_sos1
#' @export
opt_sos2 <- function(model, x, weights = NULL, name = NULL) {
  sos_add(model, x, weights, name, 2L)
}

sos_add <- function(model, x, weights, name, type) {
  check_model(model)
  sym <- resolve_symbols(model, x, "var")
  if (length(sym) < 2L) stop("A special ordered set needs at least two variables.",
                             call. = FALSE)
  weights <- weights %||% seq_along(sym)
  check_numeric(weights, "weights")
  if (length(weights) != length(sym) || anyDuplicated(weights)) {
    stop("`weights` must be distinct, one per variable.", call. = FALSE)
  }
  ord <- order(weights)
  add_special(model, "sos", list(type = type, sym = sym[ord],
                                 weights = weights[ord]), name)
}

#' @rdname opt_sos1
#' @export
opt_complementarity <- function(model, x, y, name = NULL) {
  check_model(model)
  xs <- resolve_symbols(model, x, "var")
  ys <- resolve_symbols(model, y, "var")
  if (length(xs) != length(ys)) stop("`x` and `y` must have the same length.",
                                     call. = FALSE)
  if (any(model$lower[c(xs, ys)] < 0)) {
    stop("Complementarity needs nonnegative variables (lower bounds >= 0).",
         call. = FALSE)
  }
  base <- name %||% paste0("complementarity",
                           length((model$special %||% list())$sos) + 1L)
  invisible(lapply(seq_along(xs), function(i) {
    add_special(model, "sos", list(type = 1L, sym = c(xs[i], ys[i]),
                                   weights = c(1, 2)),
                paste0(base, "[", i, "]"))
  }))
}

#' Logical constraints on binary variables
#'
#' Exact linear formulations of logical relations between binary
#' variables. R's `&&`, `||` and `!` cannot build symbolic constraints, so
#' these functions are the explicit alternative (`1 - x` is the negation of
#' a binary `x`).
#'
#' * `opt_implies(model, a, b)`: \eqn{a \Rightarrow b}, i.e. \eqn{a \le b}.
#' * `opt_and(model, z, x)`: \eqn{z = x_1 \wedge \dots \wedge x_k}.
#' * `opt_or(model, z, x)`: \eqn{z = x_1 \vee \dots \vee x_k}.
#'
#' @param model An [opt_model()].
#' @param a,b,z Binary variables (scalars, or vectors of equal length for
#'   `opt_implies()`).
#' @param x A vector of binary variables.
#' @param name Optional constraint name.
#' @return An `opt_constraint` reference, invisibly.
#' @examples
#' m <- opt_model("max")
#' x <- opt_variable(m, "x", dim = 3, type = "binary")
#' z <- opt_variable(m, "z", type = "binary")
#' opt_and(m, z, x)
#' opt_constraint(m, sum(x) <= 2)
#' opt_objective(m, z)
#' opt_objective_value(opt_solve(m))  # 0: not all three can be 1
#' @export
opt_implies <- function(model, a, b, name = NULL) {
  check_model(model)
  as_binary_symbol(model, a, "a")
  as_binary_symbol(model, b, "b")
  opt_constraint(model, a - b <= 0, name = name)
}

#' @rdname opt_implies
#' @export
opt_and <- function(model, z, x, name = NULL) {
  check_model(model)
  zs <- as_binary_symbol(model, z, "z")
  xs <- as_binary_symbol(model, x, "x")
  if (length(zs) != 1L) stop("`z` must be a single binary variable.",
                             call. = FALSE)
  k <- length(xs)
  # z <= x_i for all i, z >= sum(x) - (k - 1)
  rows <- c(rep(z, k) - x, sum(x) - z)
  opt_constraint(model, rows, lower = c(rep(-Inf, k), -Inf),
                 upper = c(rep(0, k), k - 1), name = name)
}

#' @rdname opt_implies
#' @export
opt_or <- function(model, z, x, name = NULL) {
  check_model(model)
  zs <- as_binary_symbol(model, z, "z")
  xs <- as_binary_symbol(model, x, "x")
  if (length(zs) != 1L) stop("`z` must be a single binary variable.",
                             call. = FALSE)
  k <- length(xs)
  # z >= x_i for all i, z <= sum(x)
  rows <- c(rep(z, k) - x, z - sum(x))
  opt_constraint(model, rows, lower = c(rep(0, k), -Inf),
                 upper = c(rep(Inf, k), 0), name = name)
}

# ----------------------------------------------------- standard form part

affine_matrix <- function(e, is_par, pval, col, n) {
  low <- compact_poly(lower_expr(e, is_par, pval))
  if (low$type != "poly" || length(low$qv)) {
    stop("A conic or indicator expression is not affine after parameter substitution.",
         call. = FALSE)
  }
  list(M = Matrix::sparseMatrix(i = low$li, j = col[low$lj], x = low$lv,
                                dims = c(low$len, n)),
       const = rep_len(low$const, low$len))
}

lower_special <- function(model, is_par, pval, col, n) {
  sp <- model$special %||% list()
  cones <- lapply(sp$cone %||% list(), function(cn) {
    parts <- lapply(cn$args, affine_matrix, is_par = is_par, pval = pval,
                    col = col, n = n)
    list(name = cn$name, type = cn$type, alpha = cn$alpha, k = cn$k,
         M = do.call(rbind, lapply(parts, `[[`, "M")),
         const = unlist(lapply(parts, `[[`, "const")))
  })
  ind <- list()
  for (it in sp$indicator %||% list()) {
    a <- affine_matrix(it$expr, is_par, pval, col, n)
    for (r in seq_len(nrow(a$M))) {
      row <- a$M[r, ]
      nz <- which(row != 0)
      # expression a'x + c  (sense)  0, as a'x <= -c and/or -a'x <= c
      add <- function(sign) list(name = it$name, z = col[it$z[r]],
                                 active = it$active, j = nz,
                                 v = sign * row[nz], rhs = -sign * a$const[r])
      if (it$sense %in% c("<=", "==")) ind <- c(ind, list(add(1)))
      if (it$sense %in% c(">=", "==")) ind <- c(ind, list(add(-1)))
    }
  }
  sos <- lapply(sp$sos %||% list(), function(s) {
    list(name = s$name, type = s$type, cols = col[s$sym], weights = s$weights)
  })
  list(cones = cones, indicators = ind, sos = sos)
}

special_kinds <- function(sf) {
  c(if (length(sf$cones)) "cone", if (length(sf$indicators)) "indicator",
    if (length(sf$sos)) "sos")
}

cone_types <- function(sf) unique(vapply(sf$cones, `[[`, "", "type"))

# ------------------------------------------------- exact reformulations

# Rewrite indicator and SOS constraints as big-M rows and binaries for a
# backend without native support. Columns are appended; results are
# trimmed back to the original columns and rows afterwards.
reformulate_special <- function(sf, entry) {
  native <- entry$features$special %||% character()
  need_ind <- length(sf$indicators) && !"indicator" %in% native
  need_sos <- length(sf$sos) && !"sos" %in% native
  if (!need_ind && !need_sos) return(sf)
  n0 <- length(sf$c)
  new_i <- integer(); new_j <- integer(); new_v <- numeric()
  lo <- numeric(); hi <- numeric(); nr <- 0L
  add_row <- function(j, v, l, u) {
    nr <<- nr + 1L
    new_i <<- c(new_i, rep.int(nr, length(j)))
    new_j <<- c(new_j, j)
    new_v <<- c(new_v, v)
    lo <<- c(lo, l)
    hi <<- c(hi, u)
  }
  add_col <- character()
  bound_err <- function(what) {
    stop(sprintf(paste(
      "%s needs finite bounds on its variables to be reformulated exactly",
      "for solver '%s'; add bounds or use solver = \"scip\"."), what,
      entry$name), call. = FALSE)
  }
  if (need_ind) {
    for (it in sf$indicators) {
      act <- sum(pmax(it$v * sf$lower[it$j], it$v * sf$upper[it$j]))
      if (!is.finite(act)) bound_err(sprintf("Indicator constraint '%s'", it$name))
      M <- act - it$rhs
      if (M <= 0) next # the relation always holds within the bounds
      if (it$active == 1) {
        # a'x - rhs <= M (1 - z)
        add_row(c(it$j, it$z), c(it$v, M), -Inf, it$rhs + M)
      } else {
        # a'x - rhs <= M z
        add_row(c(it$j, it$z), c(it$v, -M), -Inf, it$rhs)
      }
    }
  }
  if (need_sos) {
    for (s in sf$sos) {
      cols <- s$cols
      l <- sf$lower[cols]
      u <- sf$upper[cols]
      if (any(!is.finite(c(l, u)))) bound_err(sprintf("SOS constraint '%s'", s$name))
      k <- length(cols)
      nb <- if (s$type == 1L) k else k - 1L
      ycols <- n0 + length(add_col) + seq_len(nb)
      add_col <- c(add_col, sprintf("%s.y[%d]", s$name, seq_len(nb)))
      cover <- function(i) {
        # binaries that allow column i to be nonzero
        if (s$type == 1L) ycols[i] else ycols[intersect(c(i - 1L, i), seq_len(nb))]
      }
      for (i in seq_len(k)) {
        yc <- cover(i)
        # x_i <= u_i * sum(y)  and  x_i >= l_i * sum(y)
        add_row(c(cols[i], yc), c(1, rep(-u[i], length(yc))), -Inf, 0)
        add_row(c(cols[i], yc), c(1, rep(-l[i], length(yc))), 0, Inf)
      }
      add_row(ycols, rep(1, nb), -Inf, 1)
    }
  }
  na <- length(add_col)
  n <- n0 + na
  A0 <- sf$A
  A0 <- cbind(A0, sparse_zero(nrow(A0), na))
  A1 <- Matrix::sparseMatrix(i = new_i, j = new_j, x = new_v, dims = c(nr, n))
  out <- sf
  out$A <- rbind(A0, A1)
  out$row_lower <- c(sf$row_lower, lo)
  out$row_upper <- c(sf$row_upper, hi)
  out$row_kind <- c(sf$row_kind, rep("linear", nr))
  out$row_names <- c(sf$row_names, sprintf("reformulation[%d]", seq_len(nr)))
  out$c <- c(sf$c, numeric(na))
  if (!is.null(sf$Q)) {
    Q <- methods::as(sf$Q, "generalMatrix")
    Q <- rbind(cbind(Q, sparse_zero(n0, na)), sparse_zero(na, n))
    out$Q <- Matrix::forceSymmetric(Q, uplo = "U")
  }
  out$lower <- c(sf$lower, numeric(na))
  out$upper <- c(sf$upper, rep(1, na))
  out$type <- c(sf$type, rep("binary", na))
  out$start <- c(sf$start, rep(NA_real_, na))
  out$col_names <- c(sf$col_names, add_col)
  if (need_ind) out$indicators <- list()
  if (need_sos) out$sos <- list()
  out$reformulated <- list(columns = n0, rows = nrow(sf$A),
                           added_columns = na, added_rows = nr)
  out
}

trim_result <- function(res, sf) {
  n <- length(sf$c)
  m <- nrow(sf$A)
  if (length(res$x) > n) res$x <- res$x[seq_len(n)]
  if (length(res$col_dual) > n) res$col_dual <- res$col_dual[seq_len(n)]
  if (length(res$row_dual) > m) res$row_dual <- res$row_dual[seq_len(m)]
  res
}

# Violations of special constraints at x (columns of the standard form).
special_violations <- function(sf, x, tol) {
  out <- data.frame(constraint = character(), kind = character(),
                    violation = numeric(), stringsAsFactors = FALSE)
  add <- function(nm, kind, v) {
    out[nrow(out) + 1L, ] <<- list(nm, kind, v)
  }
  for (cn in sf$cones) {
    v <- as.vector(cn$M %*% x) + cn$const
    viol <- switch(cn$type,
      soc = sqrt(sum(v[-1L]^2)) - v[1L],
      exp = if (v[2L] > 0) v[2L] * exp(v[1L] / v[2L]) - v[3L] else
        if (v[2L] >= -tol && v[1L] <= tol && v[3L] >= -tol) 0 else
          max(-v[2L], v[1L], -v[3L]),
      pow = max(abs(v[3L]) - max(v[1L], 0)^cn$alpha * max(v[2L], 0)^(1 - cn$alpha),
                -v[1L], -v[2L]),
      psd = {
        S <- matrix(0, cn$k, cn$k)
        S[upper.tri(S, diag = TRUE)] <- v
        S[lower.tri(S)] <- t(S)[lower.tri(S)]
        -min(eigen(S, symmetric = TRUE, only.values = TRUE)$values)
      })
    add(cn$name, paste("cone", cn$type), max(viol, 0))
  }
  for (it in sf$indicators) {
    on <- abs(x[it$z] - it$active) <= 0.5
    v <- if (on) max(sum(it$v * x[it$j]) - it$rhs, 0) else 0
    add(it$name, "indicator", v)
  }
  for (s in sf$sos) {
    nzv <- abs(x[s$cols]) > max(tol, 1e-6)
    nz <- which(nzv)
    ok <- if (s$type == 1L) length(nz) <= 1L else
      length(nz) <= 2L && (length(nz) < 2L || diff(nz) == 1L)
    v <- if (ok) 0 else {
      a <- sort(abs(x[s$cols]), decreasing = TRUE)
      sum(a[-seq_len(min(length(a), s$type))])
    }
    add(s$name, paste0("sos", s$type), v)
  }
  out
}
