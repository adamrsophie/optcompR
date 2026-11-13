# Further backends: SCIP, OSQP, NLopt and ROI.

pkg_version <- function(pkg) {
  function() as.character(utils::packageVersion(pkg))
}

# Upper triangle of a one-half convention matrix as product terms:
# 0.5 x'Qx = sum_{j<k} Q_jk x_j x_k + sum_j 0.5 Q_jj x_j^2.
half_q_terms <- function(Q) {
  T <- methods::as(Matrix::triu(methods::as(Q, "generalMatrix")),
                   "TsparseMatrix")
  i <- T@i + 1L
  j <- T@j + 1L
  v <- ifelse(i == j, 0.5 * T@x, T@x)
  list(i = i, j = j, v = v)
}

row_terms <- function(A) {
  R <- methods::as(methods::as(A, "generalMatrix"), "RsparseMatrix")
  list(p = R@p, j = R@j + 1L, x = R@x)
}

# ------------------------------------------------------------------- SCIP

scip_termination <- function(s) {
  switch(s,
    optimal = "optimal", gaplimit = "optimal", infeasible = "infeasible",
    unbounded = "unbounded", infeasible_or_unbounded = "infeasible_or_unbounded",
    inforunbd = "infeasible_or_unbounded", timelimit = "time_limit",
    nodelimit = , stallnodelimit = , totalnodelimit = "node_limit",
    sollimit = , bestsollimit = "solution_limit", memlimit = "memory_limit",
    userinterrupt = "interrupted", "unknown")
}

scip_solve_sf <- function(sf, cls, control) {
  n <- length(sf$c)
  vt <- ifelse(sf$type == "binary", "B", ifelse(sf$type == "integer", "I", "C"))
  quad_obj <- !is.null(sf$Q)
  m <- scip::scip_model("optcompR")
  on.exit(try(scip::scip_model_free(m), silent = TRUE), add = TRUE)
  scip::scip_set_param(m, "display/verblevel",
                       if (isTRUE(control$verbose)) 4L else 0L)
  if (is.finite(control$time_limit)) {
    scip::scip_set_param(m, "limits/time", control$time_limit)
  }
  if (!is.null(control$mip_gap)) {
    scip::scip_set_param(m, "limits/gap", control$mip_gap)
  }
  for (nm in names(control$options)) {
    scip::scip_set_param(m, nm, control$options[[nm]])
  }
  obj <- if (quad_obj) numeric(n) else sf$c
  scip::scip_add_vars(m, obj = obj, lb = sf$lower, ub = sf$upper, vtype = vt)
  if (quad_obj) {
    # epigraph: exact reformulation, the auxiliary variable t carries the
    # objective, 0.5 x'Qx + c'x - t <= 0 (min) or >= 0 (max)
    scip::scip_add_var(m, obj = 1, lb = -Inf, ub = Inf, vtype = "C",
                       name = "optcompR_objective")
    qt <- half_q_terms(sf$Q)
    nz <- which(sf$c != 0)
    scip::scip_add_quadratic_cons(
      m, linvars = c(nz, n + 1L), lincoefs = c(sf$c[nz], -1),
      quadvars1 = qt$i, quadvars2 = qt$j, quadcoefs = qt$v,
      lhs = if (sf$sense == "max") 0 else -Inf,
      rhs = if (sf$sense == "max") Inf else 0)
  }
  rt <- row_terms(sf$A)
  for (r in seq_len(nrow(sf$A))) {
    idx <- if (rt$p[r + 1L] > rt$p[r]) (rt$p[r] + 1L):rt$p[r + 1L] else integer()
    key <- as.character(r)
    if (!is.null(sf$row_Q[[key]])) {
      qt <- half_q_terms(sf$row_Q[[key]])
      scip::scip_add_quadratic_cons(m, linvars = rt$j[idx], lincoefs = rt$x[idx],
                                    quadvars1 = qt$i, quadvars2 = qt$j,
                                    quadcoefs = qt$v, lhs = sf$row_lower[r],
                                    rhs = sf$row_upper[r])
    } else {
      scip::scip_add_linear_cons(m, vars = rt$j[idx], coefs = rt$x[idx],
                                 lhs = sf$row_lower[r], rhs = sf$row_upper[r])
    }
  }
  scip::scip_set_objective_sense(m, if (sf$sense == "max") "maximize" else
    "minimize")
  if (isTRUE(control$warm_start) && any(!is.na(sf$start))) {
    # SCIP's R binding has no API for start solutions; note it
    warn_start <- TRUE
  } else {
    warn_start <- FALSE
  }
  scip::scip_optimize(m)
  st <- scip::scip_get_status(m)
  info <- scip::scip_get_info(m)
  sol <- scip::scip_get_solution(m)
  has <- !is.null(sol$x) && length(sol$x) >= n && isTRUE(info$sol_count > 0)
  x <- if (has) sol$x[seq_len(n)] else rep(NA_real_, n)
  term <- scip_termination(st)
  obj <- if (has) {
    if (quad_obj) sum(sf$c * x) + 0.5 * sum(x * as.vector(sf$Q %*% x)) +
      sf$offset else sol$objval + sf$offset
  } else NA_real_
  list(termination = term,
       primal_status = if (has) "feasible_point" else "no_solution",
       dual_status = "not_available", x = x, objective = obj,
       row_dual = NULL, col_dual = NULL,
       iterations = c(lp = info$iterations, nodes = info$nodes),
       gap = info$gap, bound = NA_real_, message = st, local = FALSE,
       warm_start = FALSE,
       warnings = if (warn_start)
         "The scip R package has no start solution interface; start values were not used."
       else character(),
       extra = list(), raw = list(status = st, info = info, solution = sol))
}

# ------------------------------------------------------------------- OSQP

osqp_solve_sf <- function(sf, cls, control) {
  s <- if (sf$sense == "max") -1 else 1
  n <- length(sf$c)
  k <- nrow(sf$A)
  fb <- which(is.finite(sf$lower) | is.finite(sf$upper))
  B <- Matrix::sparseMatrix(i = seq_along(fb), j = fb, x = 1,
                            dims = c(length(fb), n))
  A <- rbind(methods::as(sf$A, "generalMatrix"), B)
  l <- c(sf$row_lower, sf$lower[fb])
  u <- c(sf$row_upper, sf$upper[fb])
  P <- if (is.null(sf$Q)) NULL else
    methods::as(methods::as(Matrix::triu(s * sf$Q), "generalMatrix"),
                "CsparseMatrix")
  defaults <- list(verbose = isTRUE(control$verbose), eps_abs = 1e-7,
                   eps_rel = 1e-7, polishing = TRUE, max_iter = 100000L)
  if (is.finite(control$time_limit)) defaults$time_limit <- control$time_limit
  pars <- do.call(osqp::osqpSettings, utils::modifyList(defaults,
                                                        control$options))
  res <- osqp::solve_osqp(P = P, q = s * sf$c, A = A, l = l, u = u,
                          pars = pars)
  st <- res$info$status
  term <- switch(st,
    "solved" = "optimal", "solved inaccurate" = "almost_optimal",
    "primal infeasible" = , "primal infeasible inaccurate" = "infeasible",
    # a dual infeasibility certificate does not establish primal feasibility
    "dual infeasible" = , "dual infeasible inaccurate" = "infeasible_or_unbounded",
    "maximum iterations reached" = "iteration_limit",
    "run time limit reached" = "time_limit", "unknown")
  has <- term %in% c("optimal", "almost_optimal", "iteration_limit",
                     "time_limit") && !anyNA(res$x)
  x <- if (has) res$x else rep(NA_real_, n)
  y <- -s * res$y # OSQP multipliers have the opposite sign
  rd <- if (term %in% c("optimal", "almost_optimal")) y[seq_len(k)] else NULL
  cd <- NULL
  if (!is.null(rd)) {
    cd <- numeric(n)
    cd[fb] <- y[k + seq_along(fb)]
  }
  obj <- if (has) s * res$info$obj_val + sf$offset else NA_real_
  list(termination = term,
       primal_status = if (has) "feasible_point" else "no_solution",
       dual_status = if (!is.null(rd)) "feasible_point" else "no_solution",
       x = x, objective = obj, row_dual = rd, col_dual = cd,
       iterations = c(admm = res$info$iter), gap = NA_real_,
       bound = NA_real_, message = st, local = FALSE, warm_start = FALSE,
       extra = list(), raw = res)
}

# ----------------------------------------------------------------- NLopt

nloptr_termination <- function(status) {
  if (status %in% 1:4) return("optimal")
  if (status == 5L) return("iteration_limit")
  if (status == 6L) return("time_limit")
  if (status == -4L) return("numerical_error")
  if (status == -2L) return("invalid_model")
  "unknown"
}

nloptr_solve_sf <- function(sf, cls, control) {
  s <- if (sf$sense == "max") -1 else 1
  n <- length(sf$c)
  vs <- sf$var_sym
  ns <- sf$n_sym
  obj <- sf$objective_lowered
  cons <- Filter(Negate(is.null), sf$cons_lowered)
  zfull <- function(x) {
    z <- numeric(ns)
    z[vs] <- x
    z
  }
  f <- function(x) {
    z <- zfull(x)
    jj <- jac_expr(obj, z, ns)
    list(objective = s * jj$val,
         gradient = s * as.vector(jj$J[, vs, drop = FALSE]))
  }
  # bounds on the lowered row expressions (constants included)
  lo <- sf$expr_lower
  hi <- sf$expr_upper
  eq <- is.finite(lo) & is.finite(hi) & lo == hi
  up <- is.finite(hi) & !eq
  dn <- is.finite(lo) & !eq
  gall <- function(x) {
    z <- zfull(x)
    parts <- lapply(cons, jac_expr, z = z, ns = ns)
    list(val = unlist(lapply(parts, `[[`, "val")),
         J = do.call(rbind, lapply(parts, function(p) p$J[, vs, drop = FALSE])))
  }
  g_in <- if (any(up | dn)) function(x) {
    g <- gall(x)
    J <- as.matrix(g$J)
    list(constraints = c(g$val[up] - hi[up], lo[dn] - g$val[dn]),
         jacobian = rbind(J[up, , drop = FALSE], -J[dn, , drop = FALSE]))
  } else NULL
  g_eq <- if (any(eq)) function(x) {
    g <- gall(x)
    list(constraints = g$val[eq] - lo[eq],
         jacobian = as.matrix(g$J)[eq, , drop = FALSE])
  } else NULL
  x0 <- sf$start
  miss <- is.na(x0)
  x0[miss] <- pmin(pmax(0, sf$lower[miss]), sf$upper[miss])
  algo <- control$options$algorithm %||% "NLOPT_LD_SLSQP"
  opts <- utils::modifyList(list(algorithm = algo, xtol_rel = 1e-10,
                                 ftol_rel = 1e-12, maxeval = 5000L),
                            control$options)
  if (is.finite(control$time_limit)) opts$maxtime <- control$time_limit
  if (isTRUE(control$verbose)) opts$print_level <- 1L
  res <- nloptr::nloptr(x0 = x0, eval_f = f, lb = sf$lower, ub = sf$upper,
                        eval_g_ineq = g_in, eval_g_eq = g_eq, opts = opts)
  x <- res$solution
  term <- nloptr_termination(res$status)
  primal <- "feasible_point"
  # a local solver can stop at an infeasible point: check it here
  if (length(cons)) {
    v <- gall(x)$val
    sc <- pmax(1, abs(ifelse(is.finite(lo), lo, 0)), abs(ifelse(is.finite(hi), hi, 0)))
    viol <- max(c(pmax(lo - v, v - hi, 0) / sc, 0))
    if (viol > 1e-6) {
      primal <- "infeasible_point"
      if (term == "optimal") term <- "locally_infeasible"
    }
  }
  list(termination = term, primal_status = primal,
       dual_status = "not_available", x = x,
       objective = s * res$objective + 0, row_dual = NULL, col_dual = NULL,
       iterations = c(evaluations = res$iterations), gap = NA_real_,
       bound = NA_real_, message = res$message,
       local = cls$convexity != "convex",
       warm_start = any(!miss), extra = list(), raw = res)
}

# ------------------------------------------------------------------- ROI

# ROI takes slam triplet matrices.
stm <- function(M) {
  T <- methods::as(M, "TsparseMatrix")
  slam::simple_triplet_matrix(T@i + 1L, T@j + 1L, T@x, nrow(M), ncol(M))
}

roi_solve_sf <- function(sf, cls, control) {
  n <- length(sf$c)
  objective <- if (is.null(sf$Q)) ROI::L_objective(sf$c) else
    ROI::Q_objective(Q = stm(methods::as(sf$Q, "generalMatrix")), L = sf$c)
  rows <- list()
  lo <- sf$row_lower
  hi <- sf$row_upper
  A <- methods::as(sf$A, "generalMatrix")
  eq <- is.finite(lo) & is.finite(hi) & lo == hi
  if (any(eq)) rows <- c(rows, list(ROI::L_constraint(A[eq, , drop = FALSE],
                                                      rep("==", sum(eq)),
                                                      lo[eq])))
  up <- is.finite(hi) & !eq
  if (any(up)) rows <- c(rows, list(ROI::L_constraint(stm(A[up, , drop = FALSE]),
                                                      rep("<=", sum(up)),
                                                      hi[up])))
  dn <- is.finite(lo) & !eq
  if (any(dn)) rows <- c(rows, list(ROI::L_constraint(stm(A[dn, , drop = FALSE]),
                                                      rep(">=", sum(dn)),
                                                      lo[dn])))
  con <- if (length(rows)) do.call(rbind, rows) else NULL
  bounds <- ROI::V_bound(li = seq_len(n), lb = sf$lower, ui = seq_len(n),
                         ub = sf$upper, nobj = n)
  types <- ifelse(sf$type == "binary", "B", ifelse(sf$type == "integer", "I",
                                                   "C"))
  op <- ROI::OP(objective, constraints = con, types = types, bounds = bounds,
                maximum = sf$sense == "max")
  solver <- control$options$roi_solver %||% "auto"
  ctrl <- control$options[setdiff(names(control$options), "roi_solver")]
  if (is.finite(control$time_limit)) ctrl$tm_limit <- control$time_limit * 1000
  res <- ROI::ROI_solve(op, solver = solver, control = ctrl)
  msg <- res$status$msg
  code <- res$status$code
  txt <- tolower(paste(msg$message, msg$symbol))
  term <- if (identical(code, 0L)) "optimal" else
    if (grepl("unbounded", txt) && grepl("infeas", txt)) "infeasible_or_unbounded"
    else if (grepl("infeas|nofeas", txt)) "infeasible"
    else if (grepl("unbounded|unbnd", txt)) "unbounded"
    else if (grepl("time", txt)) "time_limit" else "unknown"
  x <- if (identical(code, 0L)) ROI::solution(res, "primal") else
    rep(NA_real_, n)
  list(termination = term,
       primal_status = if (identical(code, 0L)) "feasible_point" else
         "no_solution",
       dual_status = "not_available", x = unname(x),
       objective = if (identical(code, 0L)) res$objval + sf$offset else NA_real_,
       row_dual = NULL, col_dual = NULL, iterations = NULL, gap = NA_real_,
       bound = NA_real_, message = paste(msg$solver, msg$message),
       local = FALSE, warm_start = FALSE, extra = list(), raw = res)
}

roi_version <- function() {
  as.character(utils::packageVersion("ROI"))
}

register_builtin_backends <- function() {
  register_solver(
    "scip", "scip",
    "SCIP: branch and bound with spatial branching, global for nonconvex quadratic models",
    classes = c("LP", "MILP", "QP", "MIQP", "QCQP", "MIQCQP"),
    global = c("LP", "MILP", "QP", "MIQP", "QCQP", "MIQCQP"),
    features = list(time_limit = TRUE, mip_gap = TRUE),
    solve = scip_solve_sf, version = pkg_version("scip"), priority = 20L)
  register_solver(
    "osqp", "osqp",
    "OSQP: operator splitting for convex QP and LP (first-order, moderate accuracy)",
    classes = c("LP", "QP"), convex_only = "QP", global = c("LP", "QP"),
    features = list(duals = TRUE, time_limit = TRUE, iteration_limit = TRUE),
    solve = osqp_solve_sf, version = pkg_version("osqp"), priority = 30L)
  register_solver(
    "nloptr", "nloptr",
    "NLopt SLSQP (default): local gradient-based nonlinear programming",
    classes = c("NLP", "QCQP", "QP", "LP"),
    features = list(warm_start = TRUE, time_limit = TRUE,
                    iteration_limit = TRUE),
    solve = nloptr_solve_sf, version = pkg_version("nloptr"), priority = 60L)
  register_solver(
    "roi", c("ROI", "slam"),
    "ROI: R Optimization Infrastructure with installed plugins (options$roi_solver)",
    classes = c("LP", "MILP", "QP"), convex_only = "QP",
    global = c("LP", "MILP"),
    features = list(time_limit = TRUE),
    solve = roi_solve_sf, version = roi_version, priority = 70L)
}
