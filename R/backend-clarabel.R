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

clarabel_solve_sf <- function(sf, cls, control) {
  s <- if (sf$sense == "max") -1 else 1
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
  # conic blocks in Clarabel's order: second-order, PSD, exponential, power
  by_type <- function(tp) Filter(function(cn) cn$type == tp, sf$cones)
  for (cn in by_type("soc")) {
    blocks <- c(blocks, list(-cn$M))
    rhs <- c(rhs, list(cn$const))
    cones$q <- c(cones$q, as.integer(nrow(cn$M)))
  }
  for (cn in by_type("psd")) {
    w <- svec_scale(cn$k)
    blocks <- c(blocks, list(-Matrix::Diagonal(x = w) %*% cn$M))
    rhs <- c(rhs, list(w * cn$const))
    cones$s <- c(cones$s, as.integer(cn$k))
  }
  ep <- by_type("exp")
  for (cn in ep) {
    blocks <- c(blocks, list(-cn$M))
    rhs <- c(rhs, list(cn$const))
  }
  if (length(ep)) cones$ep <- as.integer(length(ep))
  for (cn in by_type("pow")) {
    blocks <- c(blocks, list(-cn$M))
    rhs <- c(rhs, list(cn$const))
    cones$p <- c(cones$p, cn$alpha)
  }
  Ab <- methods::as(methods::as(do.call(rbind, lapply(blocks, function(b)
    methods::as(b, "generalMatrix"))), "generalMatrix"), "CsparseMatrix")
  b <- unlist(rhs)
  P <- if (is.null(sf$Q)) NULL else
    methods::as(methods::as(Matrix::triu(s * sf$Q), "generalMatrix"),
                "CsparseMatrix")
  ctrl <- utils::modifyList(list(verbose = isTRUE(control$verbose)),
                            control$options)
  if (is.finite(control$time_limit)) ctrl$time_limit <- control$time_limit
  res <- clarabel::clarabel(A = Ab, b = b, q = s * sf$c, P = P, cones = cones,
                            control = do.call(clarabel::clarabel_control, ctrl))
  term <- clarabel_termination(res$status)
  ok <- term %in% c("optimal", "almost_optimal")
  x <- if (ok) res$x else rep(NA_real_, n)
  z <- res$z
  rd <- cd <- NULL
  if (ok) {
    off <- 0L
    take <- function(k) {
      v <- z[off + seq_len(k)]
      off <<- off + k
      v
    }
    z_eq <- take(n_eq); z_u <- take(length(ur)); z_l <- take(length(lr))
    z_ub <- take(length(ub)); z_lb <- take(length(lb))
    rd <- numeric(nrow(A))
    rd[eq] <- -z_eq
    rd[ur] <- rd[ur] - z_u
    rd[lr] <- rd[lr] + z_l
    cd <- numeric(n)
    cd[ub] <- cd[ub] - z_ub
    cd[lb] <- cd[lb] + z_lb
    rd <- s * rd
    cd <- s * cd
  }
  list(termination = if (term == "almost_optimal") "almost_optimal" else term,
       primal_status = if (ok) "feasible_point" else "no_solution",
       dual_status = if (ok) "feasible_point" else "no_solution",
       x = x, objective = if (ok) s * res$obj_val + sf$offset else NA_real_,
       row_dual = rd, col_dual = cd,
       iterations = c(ipm = res$iterations), gap = NA_real_,
       bound = if (ok) s * res$obj_val_dual + sf$offset else NA_real_,
       message = names(clarabel::solver_status_descriptions())[res$status],
       local = FALSE, warm_start = FALSE,
       extra = list(cone_duals = if (ok) z[-seq_len(n_eq + n_lin)] else NULL),
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
    priority = 15L)
}
