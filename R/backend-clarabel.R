# Clarabel interior point backend for LP, convex QP and conic programs.
#
# Clarabel solves  min 0.5 x'Px + q'x  s.t.  Ax + s = b,  s in K  with dual
# z in K*. For a row a'x <= u (A row a, b = u) the optimal value moves by
# -z per unit of u; for a'x >= l (A row -a, b = -l) by +z per unit of l;
# for an equality (zero cone) by -z. Multiplying by the sign of the user's
# sense gives the shadow price d f*/d b used throughout optcompR.

clarabel_termination <- function(code) {
  switch(as.character(code),
    "2" = "optimal", "5" = "almost_optimal", "3" = "infeasible",
    "6" = "infeasible", "4" = "infeasible_or_unbounded",
    "7" = "infeasible_or_unbounded", "8" = "iteration_limit",
    "9" = "time_limit", "10" = "numerical_error", "11" = "numerical_error",
    "unknown")
}

# PSD cones use the upper triangle by columns with off-diagonal entries
# scaled by sqrt(2), so that the vector inner product equals the trace
# inner product.
svec_scale <- function(k) {
  d <- row(matrix(0, k, k)) == col(matrix(0, k, k))
  ifelse(d[upper.tri(d, diag = TRUE)], 1, sqrt(2))
}

# Conic data  A x + s = b,  s in K  shared by Clarabel and SCS: equality
# rows (zero cone), then inequality rows and finite column bounds
# (nonnegative cone), then the cones in the order soc, psd, exp, pow. SCS
# stores PSD cones by the lower triangle, Clarabel by the upper triangle.
conic_data <- function(sf, psd_lower = FALSE) {
  n <- length(sf$c)
  A <- methods::as(sf$A, "generalMatrix")
  lo <- sf$row_lower
  up <- sf$row_upper
  eq <- which(is.finite(lo) & is.finite(up) & lo == up)
  ur <- setdiff(which(is.finite(up)), eq)
  lr <- setdiff(which(is.finite(lo)), eq)
  ub <- which(is.finite(sf$upper))
  lb <- which(is.finite(sf$lower))
  I <- Matrix::Diagonal(n)
  blocks <- list(
    A[eq, , drop = FALSE], A[ur, , drop = FALSE], -A[lr, , drop = FALSE],
    I[ub, , drop = FALSE], -I[lb, , drop = FALSE])
  rhs <- list(up[eq], up[ur], -lo[lr], sf$upper[ub], -sf$lower[lb])
  n_eq <- length(eq)
  n_lin <- length(ur) + length(lr) + length(ub) + length(lb)
  cones <- list()
  if (n_eq) cones$z <- as.integer(n_eq)
  if (n_lin) cones$l <- as.integer(n_lin)
  by_type <- function(tp) Filter(function(cn) cn$type == tp, sf$cones)
  cone_rows <- list()
  add <- function(cn, M, const) {
    blocks <<- c(blocks, list(-M))
    rhs <<- c(rhs, list(const))
    cone_rows <<- c(cone_rows, list(list(name = cn$name, type = cn$type,
                                         rows = nrow(M), perm = cn$perm)))
  }
  for (cn in by_type("soc")) {
    add(cn, cn$M, cn$const)
    cones$q <- c(cones$q, as.integer(nrow(cn$M)))
  }
  for (cn in by_type("psd")) {
    w <- svec_scale(cn$k)
    perm <- if (psd_lower) psd_lower_order(cn$k) else seq_along(w)
    cn$perm <- perm
    add(cn, (Matrix::Diagonal(x = w) %*% cn$M)[perm, , drop = FALSE],
        (w * cn$const)[perm])
    cones$s <- c(cones$s, as.integer(cn$k))
  }
  ep <- by_type("exp")
  for (cn in ep) add(cn, cn$M, cn$const)
  if (length(ep)) cones$ep <- as.integer(length(ep))
  for (cn in by_type("pow")) {
    add(cn, cn$M, cn$const)
    cones$p <- c(cones$p, cn$alpha)
  }
  Ab <- methods::as(methods::as(do.call(rbind, lapply(blocks, function(b)
    methods::as(b, "generalMatrix"))), "generalMatrix"), "CsparseMatrix")
  list(A = Ab, b = unlist(rhs), cones = cones, eq = eq, ur = ur, lr = lr,
       ub = ub, lb = lb, n_eq = n_eq, n_lin = n_lin, nrow_A = nrow(A),
       cone_rows = cone_rows)
}

# Positions, in the upper-triangle-by-columns order, of the entries of the
# lower triangle by columns: entry (i, j), i >= j, equals (j, i).
psd_lower_order <- function(k) {
  up <- which(upper.tri(matrix(0, k, k), diag = TRUE))
  pos <- matrix(0L, k, k)
  pos[up] <- seq_along(up)
  lo <- which(lower.tri(matrix(0, k, k), diag = TRUE))
  t(pos)[lo]
}

# Row and column duals (shadow prices in the user's sense) from the conic
# dual vector z of  A x + s = b; and the cone part of z by cone.
conic_duals <- function(cd, z, n, s) {
  off <- 0L
  take <- function(k) {
    v <- z[off + seq_len(k)]
    off <<- off + k
    v
  }
  z_eq <- take(cd$n_eq); z_u <- take(length(cd$ur)); z_l <- take(length(cd$lr))
  z_ub <- take(length(cd$ub)); z_lb <- take(length(cd$lb))
  rd <- numeric(cd$nrow_A)
  rd[cd$eq] <- -z_eq
  rd[cd$ur] <- rd[cd$ur] - z_u
  rd[cd$lr] <- rd[cd$lr] + z_l
  col <- numeric(n)
  col[cd$ub] <- col[cd$ub] - z_ub
  col[cd$lb] <- col[cd$lb] + z_lb
  cone <- list()
  for (cr in cd$cone_rows) {
    v <- take(cr$rows)
    # back to the upper-triangle order used by optcompR
    if (!is.null(cr$perm)) v[cr$perm] <- v
    cone[[cr$name]] <- v
  }
  list(row = s * rd, col = s * col, cone = cone,
       cone_vector = z[seq_along(z) > cd$n_eq + cd$n_lin])
}

clarabel_solve_sf <- function(sf, cls, control) {
  s <- if (sf$sense == "max") -1 else 1
  n <- length(sf$c)
  cd <- conic_data(sf)
  Ab <- cd$A
  b <- cd$b
  cones <- cd$cones
  P <- if (is.null(sf$Q)) NULL else
    methods::as(Matrix::forceSymmetric(s * sf$Q, uplo = "U"), "CsparseMatrix")
  ctrl <- utils::modifyList(list(verbose = isTRUE(control$verbose)),
                            control$options)
  if (is.finite(control$time_limit)) ctrl$time_limit <- control$time_limit
  res <- clarabel::clarabel(A = Ab, b = b, q = s * sf$c, P = P, cones = cones,
                            control = do.call(clarabel::clarabel_control, ctrl))
  term <- clarabel_termination(res$status)
  ok <- term %in% c("optimal", "almost_optimal")
  x <- if (ok) res$x else rep(NA_real_, n)
  du <- if (ok) conic_duals(cd, res$z, n, s) else NULL
  list(termination = if (term == "almost_optimal") "almost_optimal" else term,
       primal_status = if (ok) "feasible_point" else "no_solution",
       dual_status = if (ok) "feasible_point" else "no_solution",
       x = x, objective = if (ok) s * res$obj_val + sf$offset else NA_real_,
       row_dual = du$row, col_dual = du$col,
       iterations = c(ipm = res$iterations), gap = NA_real_,
       bound = if (ok) s * res$obj_val_dual + sf$offset else NA_real_,
       message = names(clarabel::solver_status_descriptions())[res$status],
       local = FALSE, warm_start = FALSE,
       extra = list(cone_duals = du$cone_vector, cone_dual = du$cone),
       raw = res)
}

register_clarabel <- function() {
  register_solver(
    "clarabel", "clarabel",
    "Clarabel: interior point for LP, convex QP and SOC, PSD, exponential and power cone programs",
    classes = c("LP", "QP", "SOCP", "CONIC"), convex_only = c("QP", "SOCP", "CONIC"),
    global = c("LP", "QP", "SOCP", "CONIC"),
    features = list(duals = TRUE, time_limit = TRUE, iteration_limit = TRUE,
                    cones = c("soc", "psd", "exp", "pow")),
    solve = clarabel_solve_sf, version = pkg_version("clarabel"),
    priority = 15L,
    # its interior point method is much faster than the HiGHS active set QP
    # solver beyond a few hundred variables (0.11 s against 0.44 s at 1,000
    # and 6 s against 110 s at 5,000 in benchmarks/run.R)
    first_for = "QP")
}
