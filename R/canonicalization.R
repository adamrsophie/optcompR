# Lowering a model to a numerical standard form.
#
#   minimize or maximize  0.5 x'Qx + c'x + offset   (or a nonlinear f(x))
#   subject to            row_lower <= A x + q_i(x) <= row_upper
#                         lower <= x <= upper,  x_j integer for j in J
#
# Parameters are substituted first. Q uses the one-half convention, so the
# quadratic term v * x_j * x_k (j != k) gives Q[j, k] = Q[k, j] = v and
# v * x_j^2 gives Q[j, j] = 2 v. Constant terms of rows move to the bounds;
# the objective constant is kept in `offset`.

sparse_zero <- function(nr, nc) {
  Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(),
                       dims = c(nr, nc))
}

half_q <- function(j, k, v, n) {
  Matrix::forceSymmetric(
    Matrix::sparseMatrix(i = c(j, k), j = c(k, j), x = c(v, v),
                         dims = c(n, n)), uplo = "U")
}

#' Standard form of a model
#'
#' Substitutes parameter values and lowers the model to sparse matrices in
#' the form \eqn{\min / \max\ \frac12 x^\top Q x + c^\top x + \mathrm{offset}}
#' subject to \eqn{l_r \le A x + q_i(x) \le u_r}, \eqn{l \le x \le u} and
#' integrality. Rows keep a mapping back to the user's constraint blocks.
#' Results are cached and the cache is invalidated by any change to the model
#' structure, bounds or parameter values.
#'
#' @param model An [opt_model()].
#' @return An `opt_standard_form` list with elements `sense`, `c`, `Q`,
#'   `offset`, `A`, `row_lower`, `row_upper`, `lower`, `upper`, `type`,
#'   `col_names`, `row_names`, `row_block`, `row_kind` (`"linear"`,
#'   `"quadratic"` or `"nonlinear"`), `row_Q` (list of row Hessian matrices
#'   in the one-half convention for quadratic rows), and the lowered
#'   expressions used by nonlinear solvers.
#' @examples
#' m <- opt_model("max")
#' x <- opt_variable(m, "x", lower = 0)
#' y <- opt_variable(m, "y", lower = 0)
#' opt_constraint(m, x + y <= 4)
#' opt_objective(m, 3 * x + 2 * y + 1)
#' sf <- opt_standard_form(m)
#' sf$c
#' as.matrix(sf$A)
#' sf$offset
#' @export
opt_standard_form <- function(model) {
  check_model(model)
  key <- paste(model$structure_version, model$parameter_version, sep = ":")
  if (identical(model$cache$sf_key, key)) return(model$cache$sf)
  sf <- build_standard_form(model)
  model$cache$sf <- sf
  model$cache$sf_key <- key
  sf
}

build_standard_form <- function(model) {
  ns <- n_symbols(model)
  is_par <- model$sym_kind == "par"
  pval <- model$value
  vs <- which(!is_par)
  n <- length(vs)
  col <- rep(NA_integer_, ns)
  col[vs] <- seq_len(n)

  # objective
  obj <- model$objective
  obj_low <- if (is.null(obj)) as_expr(0) else lower_expr(obj, is_par, pval)
  cvec <- numeric(n)
  Q <- NULL
  offset <- 0
  obj_kind <- "linear"
  if (obj_low$type == "poly") {
    o <- compact_poly(obj_low)
    cvec <- acc_sum(col[o$lj], o$lv, n)
    offset <- o$const
    if (length(o$qv)) {
      Q <- half_q(col[o$qj], col[o$qk], o$qv, n)
      obj_kind <- "quadratic"
    }
  } else {
    obj_kind <- "nonlinear"
    offset <- NA_real_
  }

  # constraints: collect per block, concatenate once
  nb <- length(model$cons)
  parts <- vector("list", nb)
  row_Q <- list()
  cons_low <- vector("list", nb)
  off <- 0L
  for (b in seq_len(nb)) {
    blk <- model$cons[[b]]
    e <- lower_expr(blk$expr, is_par, pval)
    cons_low[[b]] <- e
    m <- e$len
    kind <- rep.int("linear", m)
    if (e$type == "poly") {
      if (length(e$qv)) e <- compact_poly(e)
      parts[[b]] <- list(i = off + e$li, j = col[e$lj], v = e$lv,
                         lo = blk$lower - e$const, hi = blk$upper - e$const)
      if (length(e$qv)) {
        qrows <- sort(unique(e$qi))
        kind[qrows] <- "quadratic"
        for (r in qrows) {
          t <- e$qi == r
          row_Q[[as.character(off + r)]] <- half_q(col[e$qj[t]], col[e$qk[t]],
                                                  e$qv[t], n)
        }
      }
    } else {
      kind[] <- "nonlinear"
      parts[[b]] <- list(i = integer(), j = integer(), v = numeric(),
                         lo = blk$lower, hi = blk$upper)
    }
    parts[[b]]$kind <- kind
    parts[[b]]$block <- rep.int(b, m)
    parts[[b]]$elem <- seq_len(m)
    off <- off + m
  }
  pull <- function(f, empty) {
    if (nb == 0L) return(empty)
    unlist(lapply(parts, `[[`, f), use.names = FALSE) %||% empty
  }
  ri <- pull("i", integer()); rj <- pull("j", integer()); rv <- pull("v", numeric())
  rlo <- pull("lo", numeric()); rhi <- pull("hi", numeric())
  rkind <- pull("kind", character()); rblock <- pull("block", integer())
  relem <- pull("elem", integer())
  A <- Matrix::sparseMatrix(i = ri, j = rj, x = rv, dims = c(off, n))
  row_names <- unlist(lapply(model$cons, function(blk) {
    element_labels(blk$name, blk$expr$len, blk$expr$dim)
  }), use.names = FALSE) %||% character()

  structure(list(
    sense = model$sense, c = cvec, Q = Q, offset = offset,
    objective_kind = obj_kind,
    A = A, row_lower = rlo, row_upper = rhi, row_block = rblock,
    row_element = relem, row_kind = rkind, row_Q = row_Q,
    row_names = row_names,
    lower = model$lower[vs], upper = model$upper[vs], type = model$type[vs],
    start = model$start[vs],
    col_names = model$sym_name[vs], var_sym = vs, col_of_sym = col,
    n_sym = ns, objective_lowered = obj_low, cons_lowered = cons_low,
    parameter_values = pval,
    key = paste(model$structure_version, model$parameter_version, sep = ":")
  ), class = "opt_standard_form")
}

#' @export
print.opt_standard_form <- function(x, ...) {
  cat("<opt_standard_form>\n")
  cat(sprintf("  %s, %d columns, %d rows, %d nonzeros\n", x$sense,
              length(x$c), nrow(x$A), length(x$A@x)))
  cat(sprintf("  objective: %s\n", x$objective_kind))
  tab <- table(x$row_kind)
  if (length(tab)) cat("  rows:", paste(names(tab), tab, sep = " = ",
                                        collapse = ", "), "\n")
  invisible(x)
}
