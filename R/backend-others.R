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
    optimal = "optimal", gaplimit = "gap_limit", infeasible = "infeasible",
    unbounded = "unbounded", infeasible_or_unbounded = "infeasible_or_unbounded",
    inforunbd = "infeasible_or_unbounded", timelimit = "time_limit",
    nodelimit = , stallnodelimit = , totalnodelimit = "node_limit",
    sollimit = , bestsollimit = "solution_limit", memlimit = "memory_limit",
    userinterrupt = "interrupted", "unknown")
}

# SCIP reports the gap |p - d| / min(|p|, |d|) but not the dual bound d;
# invert it when p and d have the same sign (otherwise the gap is infinite).
scip_dual_bound <- function(p, gap, sense) {
  if (!is.finite(p) || !is.finite(gap %||% NA_real_) || gap >= 1e19) {
    return(NA_real_)
  }
  if (gap == 0) return(p)
  # minimization: d <= p; maximization: d >= p
  toward_zero <- (sense == "min") == (p > 0)
  if (toward_zero) p / (1 + gap) else p * (1 + gap)
}

scip_solve_sf <- function(sf, cls, control) {
  n_model <- length(sf$c)
  if (cls$class %in% c("NLP", "MINLP")) {
    why <- lift_refusal(sf)
    if (!is.null(why)) stop(why, ".", call. = FALSE)
    sf <- lift_sf(sf)
  }
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
  qg <- row_quad_groups(sf)
  rq <- sf$row_quad
  for (r in seq_len(nrow(sf$A))) {
    idx <- if (rt$p[r + 1L] > rt$p[r]) (rt$p[r] + 1L):rt$p[r + 1L] else integer()
    qi <- qg[[as.character(r)]]
    if (length(qi)) {
      scip::scip_add_quadratic_cons(m, linvars = rt$j[idx], lincoefs = rt$x[idx],
                                    quadvars1 = rq$j[qi], quadvars2 = rq$k[qi],
                                    quadcoefs = rq$v[qi], lhs = sf$row_lower[r],
                                    rhs = sf$row_upper[r])
    } else {
      scip::scip_add_linear_cons(m, vars = rt$j[idx], coefs = rt$x[idx],
                                 lhs = sf$row_lower[r], rhs = sf$row_upper[r])
    }
  }
  scip_add_special(m, sf, n + quad_obj)
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
  # lifted auxiliary variables are not model variables
  x <- x[seq_len(n_model)]
  list(termination = term,
       primal_status = if (has) "feasible_point" else "no_solution",
       dual_status = "not_available", x = x, objective = obj,
       row_dual = NULL, col_dual = NULL,
       iterations = c(lp = info$iterations, nodes = info$nodes),
       gap = info$gap,
       bound = if (has) scip_dual_bound(obj - sf$offset, info$gap, sf$sense) +
         sf$offset else NA_real_,
       message = st, local = FALSE,
       warm_start = FALSE,
       warnings = if (warn_start)
         "The scip R package has no start solution interface; start values were not used."
       else character(),
       extra = list(), raw = list(status = st, info = info, solution = sol))
}

# SCIP takes nonlinear models only when they lift exactly to quadratic
# rows. In automatic selection, convex nonlinear models go to the local NLP
# solvers, whose optima are then global.
scip_check <- function(sf, cls, auto) {
  if (!cls$class %in% c("NLP", "MINLP")) return(NULL)
  if (auto && cls$convexity == "convex") {
    return("scip is not chosen automatically for convex nonlinear models")
  }
  lift_refusal(sf)
}

# Cones, indicators and SOS constraints for SCIP. `nvars` counts the
# variables already in the SCIP model.
scip_add_special <- function(m, sf, nvars) {
  for (cn in sf$cones) {
    # ||X x + b|| <= t'x + c  as  ||X x + b||^2 - (t'x + c)^2 <= 0 with
    # t'x + c >= 0, which is the same set
    M <- methods::as(cn$M, "generalMatrix")
    tr <- M[1L, , drop = FALSE]
    X <- M[-1L, , drop = FALSE]
    b <- cn$const[-1L]
    c0 <- cn$const[1L]
    H <- Matrix::crossprod(X) - Matrix::crossprod(tr)
    lin <- 2 * as.vector(Matrix::crossprod(X, b)) - 2 * c0 * as.vector(tr)
    T <- methods::as(Matrix::triu(methods::as(H, "generalMatrix")),
                     "TsparseMatrix")
    i <- T@i + 1L
    j <- T@j + 1L
    v <- ifelse(i == j, T@x, 2 * T@x)
    keep <- v != 0
    nz <- which(lin != 0)
    scip::scip_add_quadratic_cons(m, linvars = nz, lincoefs = lin[nz],
                                  quadvars1 = i[keep], quadvars2 = j[keep],
                                  quadcoefs = v[keep], lhs = -Inf,
                                  rhs = c0^2 - sum(b^2))
    tn <- which(as.vector(tr) != 0)
    scip::scip_add_linear_cons(m, vars = tn, coefs = as.vector(tr)[tn],
                               lhs = -c0, rhs = Inf)
  }
  for (it in sf$indicators) {
    z <- it$z
    if (it$active == 0) {
      # indicator on the complement w = 1 - z
      scip::scip_add_var(m, obj = 0, lb = 0, ub = 1, vtype = "B")
      nvars <- nvars + 1L
      scip::scip_add_linear_cons(m, vars = c(z, nvars), coefs = c(1, 1),
                                 lhs = 1, rhs = 1)
      z <- nvars
    }
    scip::scip_add_indicator_cons(m, binvar = z, vars = it$j, coefs = it$v,
                                  rhs = it$rhs)
  }
  for (s in sf$sos) {
    if (s$type == 1L) scip::scip_add_sos1_cons(m, vars = s$cols, weights = s$weights)
    else scip::scip_add_sos2_cons(m, vars = s$cols, weights = s$weights)
  }
  invisible(nvars)
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
  settings <- utils::modifyList(defaults, control$options)
  pars <- do.call(osqp::osqpSettings, settings)
  # the OSQP workspace is kept per model and updated in place (costs,
  # bounds and matrix values) while the structure, the sparsity patterns
  # and the settings stay the same; OSQP then warm starts from its last
  # solution and keeps its factorization when only q, l and u change
  A <- methods::as(A, "CsparseMatrix")
  cache <- control$cache
  data <- list(q = s * sf$c, l = l, u = u, Ax = A@x,
               Px = if (is.null(P)) numeric() else P@x)
  sig <- list(key = sf$structure_key, n = n, k = k, fb = fb, Ai = A@i,
              Ap = A@p, Pi = if (is.null(P)) NULL else P@i,
              Pp = if (is.null(P)) NULL else P@p, settings = settings)
  sess <- if (isTRUE(control$warm_start) && is.environment(cache)) cache$osqp
  reused <- !is.null(sess) && identical(sess$sig, sig) &&
    osqp_update(sess, data)
  work <- if (reused) sess$work else
    osqp::osqp(P = P, q = data$q, A = A, l = l, u = u, pars = pars)
  if (is.environment(cache)) {
    cache$osqp <- list(work = work, sig = sig, data = data)
  }
  res <- osqp_method(work, "Solve")()
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
  # an iterate at a limit is a solution only if it is feasible
  feas <- has && term %in% c("optimal", "almost_optimal") ||
    has && sf_point_feasible(sf, x)
  y <- -s * res$y # OSQP multipliers have the opposite sign
  rd <- if (term %in% c("optimal", "almost_optimal")) y[seq_len(k)] else NULL
  cd <- NULL
  if (!is.null(rd)) {
    cd <- numeric(n)
    cd[fb] <- y[k + seq_along(fb)]
  }
  obj <- if (has) s * res$info$obj_val + sf$offset else NA_real_
  list(termination = term,
       primal_status = if (feas) "feasible_point" else if (has)
         "infeasible_point" else "no_solution",
       dual_status = if (!is.null(rd)) "feasible_point" else "no_solution",
       x = x, objective = obj, row_dual = rd, col_dual = cd,
       iterations = c(admm = res$info$iter), gap = NA_real_,
       bound = NA_real_, message = st, local = FALSE, warm_start = reused,
       extra = list(reused = reused), raw = res)
}

# A method of an OSQP workspace: a property of the S7 object in osqp 1.0,
# an element of the R6 object before.
osqp_method <- function(work, name) {
  if (inherits(work, "S7_object")) {
    getExportedValue("S7", "prop")(work, name)
  } else {
    work[[name]]
  }
}

# Pass the changed data to a kept OSQP workspace; FALSE when it cannot be
# used any more, so a new one is set up.
osqp_update <- function(sess, data) {
  old <- sess$data
  args <- list()
  if (!identical(old$q, data$q)) args$q <- data$q
  if (!identical(old$l, data$l)) args$l <- data$l
  if (!identical(old$u, data$u)) args$u <- data$u
  if (!identical(old$Ax, data$Ax)) args$Ax <- data$Ax
  if (!identical(old$Px, data$Px)) args$Px <- data$Px
  # a workspace restored from a saved model has a dead native pointer;
  # probe it even when nothing changed
  live <- tryCatch({
    osqp_method(sess$work, "GetDims")()
    TRUE
  }, error = function(e) FALSE)
  if (!live) return(FALSE)
  if (!length(args)) return(TRUE)
  tryCatch({
    do.call(osqp_method(sess$work, "Update"), args)
    TRUE
  }, error = function(e) FALSE)
}

# ----------------------------------------------------------------- NLopt

nloptr_termination <- function(status) {
  # 2 is NLOPT_STOPVAL_REACHED: the objective target was met, which says
  # nothing about optimality
  if (status == 2L) return("objective_limit")
  if (status %in% c(1L, 3L, 4L)) return("optimal")
  if (status == 5L) return("iteration_limit")
  if (status == 6L) return("time_limit")
  if (status == -4L) return("numerical_error")
  if (status == -2L) return("invalid_model")
  "unknown"
}

nloptr_solve_sf <- function(sf, cls, control) {
  s <- if (sf$sense == "max") -1 else 1
  n <- length(sf$c)
  # SLSQP can cycle without end on an infeasible problem; with linear rows
  # an LP proves infeasibility first
  if (linear_infeasible(sf)) {
    return(list(termination = "infeasible", primal_status = "no_solution",
                dual_status = "not_available", x = rep(NA_real_, n),
                objective = NA_real_, row_dual = NULL, col_dual = NULL,
                iterations = NULL, gap = NA_real_, bound = NA_real_,
                message = "linear rows and bounds are infeasible (HiGHS LP)",
                local = FALSE, warm_start = FALSE, extra = list(), raw = NULL))
  }
  vs <- sf$var_sym
  n0 <- length(vs)
  ns <- sf$n_sym
  obj <- sf$objective_lowered
  cons <- Filter(Negate(is.null), sf$cons_lowered)
  # columns and linear rows added by reformulate_special()
  extra <- if (is.null(sf$reformulated)) integer() else
    seq.int(sf$reformulated$rows + 1L, length.out = sf$reformulated$added_rows)
  E <- as.matrix(sf$A[extra, , drop = FALSE])
  pad <- function(J) if (n > n0) cbind(J, matrix(0, nrow(J), n - n0)) else J
  zfull <- function(x) {
    z <- numeric(ns)
    z[vs] <- x[seq_len(n0)]
    z
  }
  algo <- control$options$algorithm %||% "NLOPT_LD_SLSQP"
  # derivative-free algorithms (NLOPT_LN_*, NLOPT_GN_*) get values only
  dfree <- grepl("^NLOPT_[GL]N_", algo)
  f <- if (dfree) function(x) s * eval_expr(obj, zfull(x)) else function(x) {
    z <- zfull(x)
    jj <- jac_expr(obj, z, ns)
    list(objective = s * jj$val,
         gradient = c(s * as.vector(jj$J[, vs, drop = FALSE]),
                      numeric(n - n0)))
  }
  # bounds on the lowered row expressions (constants included)
  lo <- c(sf$expr_lower, sf$row_lower[extra])
  hi <- c(sf$expr_upper, sf$row_upper[extra])
  eq <- is.finite(lo) & is.finite(hi) & lo == hi
  up <- is.finite(hi) & !eq
  dn <- is.finite(lo) & !eq
  # NLopt asks for the inequality and equality parts at the same point:
  # the rows are evaluated once per point
  last_x <- NULL
  last_g <- NULL
  gall <- function(x) {
    if (identical(x, last_x)) return(last_g)
    last_x <<- x
    last_g <<- gall_at(x)
    last_g
  }
  gall_at <- function(x) {
    z <- zfull(x)
    if (dfree) {
      return(list(val = c(unlist(lapply(cons, eval_expr, z = z)),
                          as.vector(E %*% x)), J = NULL))
    }
    parts <- lapply(cons, jac_expr, z = z, ns = ns)
    J <- if (length(parts)) pad(as.matrix(do.call(rbind, lapply(parts,
      function(p) p$J[, vs, drop = FALSE])))) else matrix(0, 0, n)
    list(val = c(unlist(lapply(parts, `[[`, "val")), as.vector(E %*% x)),
         J = rbind(J, E))
  }
  g_in <- if (any(up | dn)) function(x) {
    g <- gall(x)
    if (dfree) return(c(g$val[up] - hi[up], lo[dn] - g$val[dn]))
    J <- as.matrix(g$J)
    list(constraints = c(g$val[up] - hi[up], lo[dn] - g$val[dn]),
         jacobian = rbind(J[up, , drop = FALSE], -J[dn, , drop = FALSE]))
  } else NULL
  g_eq <- if (any(eq)) function(x) {
    g <- gall(x)
    if (dfree) return(g$val[eq] - lo[eq])
    list(constraints = g$val[eq] - lo[eq],
         jacobian = as.matrix(g$J)[eq, , drop = FALSE])
  } else NULL
  st <- default_start(sf, control)
  x0 <- st$x
  miss <- st$miss
  ok0 <- tryCatch({
    f0 <- f(x0)
    all(is.finite(if (dfree) f0 else c(f0$objective, f0$gradient))) &&
      (!(length(cons) || length(extra)) || all(is.finite(gall(x0)$val)))
  }, error = function(e) FALSE)
  if (!ok0) stop_nonfinite_start()
  opts <- utils::modifyList(list(algorithm = algo, xtol_rel = 1e-10,
                                 ftol_rel = 1e-12, maxeval = 5000L),
                            control$options)
  if (is.finite(control$time_limit)) opts$maxtime <- control$time_limit
  if (isTRUE(control$verbose)) opts$print_level <- 1L
  res <- nloptr::nloptr(x0 = x0, eval_f = f, lb = sf$lower, ub = sf$upper,
                        eval_g_ineq = g_in, eval_g_eq = g_eq, opts = opts)
  x <- res$solution
  term <- nloptr_termination(res$status)
  primal <- if (anyNA(x) || any(!is.finite(x))) "no_solution" else
    "feasible_point"
  # a local solver can stop at an infeasible point: check it here
  if (primal == "feasible_point" && (length(cons) || length(extra))) {
    v <- gall(x)$val
    viol <- max(c(bound_violation(v, lo, hi), 0))
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
  if (any(eq)) rows <- c(rows, list(ROI::L_constraint(stm(A[eq, , drop = FALSE]),
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
  # ROI registers its plugins when attached; load them explicitly here
  plugins <- if (identical(solver, "auto")) {
    unname(ROI::ROI_installed_solvers())
  } else {
    paste0("ROI.plugin.", solver)
  }
  for (p in plugins) requireNamespace(p, quietly = TRUE)
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
    "SCIP: branch and bound with spatial branching, global for nonconvex quadratic and polynomial models",
    classes = c("LP", "MILP", "QP", "MIQP", "QCQP", "MIQCQP", "SOCP", "MISOCP",
                "NLP", "MINLP"),
    global = c("LP", "MILP", "QP", "MIQP", "QCQP", "MIQCQP", "SOCP", "MISOCP",
               "NLP", "MINLP"),
    features = list(time_limit = TRUE, mip_gap = TRUE, cones = "soc",
                    special = c("indicator", "sos")),
    solve = scip_solve_sf, version = pkg_version("scip"), priority = 20L,
    check = scip_check)
  register_solver(
    "osqp", "osqp",
    "OSQP: operator splitting for convex QP and LP (first-order, moderate accuracy)",
    classes = c("LP", "QP"), convex_only = "QP", global = c("LP", "QP"),
    features = list(duals = TRUE, time_limit = TRUE, iteration_limit = TRUE),
    solve = osqp_solve_sf, version = pkg_version("osqp"), priority = 30L,
    feasibility_tol = 1e-2, accuracy_options = c("eps_abs", "eps_rel"))
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
