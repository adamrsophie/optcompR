# SCS operator splitting backend for LP, convex QP and conic programs.
#
# SCS solves  min 0.5 x'Px + c'x  s.t.  Ax + s = b,  s in K  with dual y in
# K*, the same form as Clarabel, so the conic data and the dual mapping are
# shared (conic_data(), conic_duals()). SCS is a first-order method: its
# solutions are accurate to about 1e-4 relative by default.

# A dual infeasibility certificate (code -1) shows that the problem is
# unbounded only when it is primal feasible, which SCS does not establish,
# so it is reported as infeasible_or_unbounded. An inaccurate solution
# (code 2) stopped by max_iters or the time limit is reported as that limit.
scs_termination <- function(code, status = "") {
  if (identical(as.character(code), "2")) {
    if (grepl("max_iters", status, fixed = TRUE)) return("iteration_limit")
    if (grepl("time", status, fixed = TRUE)) return("time_limit")
  }
  switch(as.character(code),
    "1" = "optimal", "2" = "almost_optimal",
    "-1" = "infeasible_or_unbounded",
    "-2" = "infeasible", "-6" = "unknown", "-7" = "unknown",
    "-4" = "numerical_error", "-5" = "interrupted", "unknown")
}

scs_solve_sf <- function(sf, cls, control) {
  s <- if (sf$sense == "max") -1 else 1
  n <- length(sf$c)
  cd <- conic_data(sf, psd_lower = TRUE)
  P <- if (is.null(sf$Q)) NULL else
    methods::as(Matrix::forceSymmetric(s * sf$Q, uplo = "U"), "CsparseMatrix")
  ctrl <- utils::modifyList(list(verbose = isTRUE(control$verbose)),
                            control$options)
  if (is.finite(control$time_limit)) {
    ctrl$time_limit_secs <- control$time_limit
  }
  initial <- NULL
  if (isTRUE(control$warm_start) && !anyNA(sf$start)) {
    initial <- list(x = sf$start, y = numeric(nrow(cd$A)),
                    s = as.vector(cd$b - cd$A %*% sf$start))
    ctrl$warm_start <- TRUE
  }
  res <- scs::scs(A = cd$A, b = cd$b, obj = s * sf$c, P = P, cone = cd$cones,
                  initial = initial,
                  control = do.call(scs::scs_control, ctrl))
  code <- as.character(res$info$status_val)
  term <- scs_termination(code, res$info$status)
  # an inaccurate iterate (code 2) is a point to check, not a solution:
  # opt_solve() checks it against the model, and its dual objective is no
  # bound
  exact <- code == "1"
  ok <- code %in% c("1", "2")
  x <- if (ok) res$x else rep(NA_real_, n)
  du <- if (ok) conic_duals(cd, res$y, n, s) else NULL
  list(termination = term,
       primal_status = if (exact) "feasible_point" else if (ok) "unknown" else
         "no_solution",
       dual_status = if (exact) "feasible_point" else if (ok) "unknown" else
         "no_solution",
       x = x, objective = if (ok) s * res$info$pobj + sf$offset else NA_real_,
       row_dual = du$row, col_dual = du$col,
       iterations = c(admm = res$info$iter), gap = NA_real_,
       bound = if (exact) s * res$info$dobj + sf$offset else NA_real_,
       message = res$info$status, local = FALSE,
       warm_start = !is.null(initial),
       extra = list(cone_duals = du$cone_vector, cone_dual = du$cone),
       raw = res)
}

register_scs <- function() {
  register_solver(
    "scs", "scs",
    "SCS: operator splitting for LP, convex QP and SOC, PSD, exponential and power cone programs (first-order, moderate accuracy)",
    classes = c("LP", "QP", "SOCP", "CONIC"),
    convex_only = c("QP", "SOCP", "CONIC"),
    global = c("LP", "QP", "SOCP", "CONIC"),
    features = list(duals = TRUE, time_limit = TRUE, iteration_limit = TRUE,
                    warm_start = TRUE, cones = c("soc", "psd", "exp", "pow")),
    solve = scs_solve_sf, version = pkg_version("scs"), priority = 35L,
    feasibility_tol = 1e-2)
}
