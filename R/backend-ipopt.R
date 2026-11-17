# IPOPT interior point backend through the optional C binding in
# src/ipopt.c. The Jacobian and Hessian structures passed to IPOPT are the
# sparsity patterns of the expression graph.

#' Is the IPOPT backend available?
#'
#' optcompR compiles an optional binding to the IPOPT C interface when
#' `pkg-config` finds IPOPT at installation (for example after
#' `brew install ipopt` or `apt install coinor-libipopt-dev`). Reinstall
#' optcompR after installing IPOPT to enable `solver = "ipopt"`.
#'
#' @return `TRUE` when the IPOPT backend was compiled in.
#' @examples
#' opt_ipopt_available()
#' @export
opt_ipopt_available <- function() {
  isTRUE(.Call(optcompR_ipopt_available) > 0L)
}

ipopt_termination <- function(status) {
  switch(as.character(status),
    "0" = "optimal", "1" = "almost_optimal",
    "2" = "locally_infeasible",
    # diverging iterates also occur on bounded problems
    "4" = "unknown",
    "-1" = "iteration_limit", "-4" = "time_limit", "-5" = "time_limit",
    "3" = , "-2" = , "-3" = , "-13" = "numerical_error",
    "-10" = , "-11" = , "-12" = "invalid_model",
    "unknown")
}

ipopt_solve_sf <- function(sf, cls, control) {
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
  # IPOPT asks for values and first derivatives at the same point in
  # separate callbacks; each is computed once per point
  memo <- function(fun) {
    lx <- NULL
    lv <- NULL
    function(x) {
      if (!identical(x, lx)) {
        lv <<- fun(x)
        lx <<- x
      }
      lv
    }
  }
  objp <- memo(function(x) jac_expr(obj, zfull(x), ns))
  f <- function(x) s * objp(x)$val
  grad <- function(x) s * as.vector(objp(x)$J[, vs, drop = FALSE])
  parts <- memo(function(x) lapply(cons, jac_expr, z = zfull(x), ns = ns))
  g <- function(x) as.numeric(unlist(lapply(parts(x), `[[`, "val")))
  lens <- vapply(cons, function(e) e$len, 0L)
  # sparse structures: structural nonzeros of the Jacobian rows and of the
  # lower triangle of the Lagrangian Hessian, in the variable columns
  m <- sum(lens)
  JP <- if (m) do.call(rbind, lapply(cons, function(e)
    jac_pattern(e, ns)[, vs, drop = FALSE])) else sparse_zero(0L, n)
  JP <- methods::as(JP, "TsparseMatrix")
  jidx <- cbind(JP@i + 1L, JP@j + 1L)
  HP <- hess_pattern(obj, ns)
  for (e in cons) HP <- HP + hess_pattern(e, ns)
  HP <- methods::as(Matrix::tril(HP[vs, vs, drop = FALSE]), "TsparseMatrix")
  hidx <- cbind(HP@i + 1L, HP@j + 1L)
  jac <- function(x) {
    J <- do.call(rbind, lapply(parts(x), function(p) p$J[, vs, drop = FALSE]))
    as.numeric(J[jidx])
  }
  hess <- function(x, obj_factor, lambda) {
    z <- zfull(x)
    H <- s * obj_factor * hess_expr(obj, z, ns, 1)
    off <- 0L
    for (k in seq_along(cons)) {
      w <- lambda[off + seq_len(lens[k])]
      if (any(w != 0)) H <- H + hess_expr(cons[[k]], z, ns, w)
      off <- off + lens[k]
    }
    as.numeric(H[vs, vs, drop = FALSE][hidx])
  }
  st <- default_start(sf, control)
  x0 <- st$x
  miss <- st$miss
  ok0 <- tryCatch(all(is.finite(c(f(x0), grad(x0), g(x0)))),
                  error = function(e) FALSE)
  if (!ok0) stop_nonfinite_start()
  # external functions without a Hessian: quasi-Newton approximation
  exact_h <- tryCatch({
    hess(x0, 1, rep(1, sum(lens)))
    TRUE
  }, error = function(e) FALSE)
  big <- function(v) pmin(pmax(v, -1e300), 1e300)
  opts <- list(print_level = if (isTRUE(control$verbose)) 5L else 0L,
               sb = "yes", tol = 1e-9)
  if (is.finite(control$time_limit)) {
    # max_wall_time exists from IPOPT 3.14; older versions have CPU time only
    if (.Call(optcompR_ipopt_available) >= 2L) {
      opts$max_wall_time <- control$time_limit
    } else {
      opts$max_cpu_time <- control$time_limit
    }
  }
  opts <- utils::modifyList(opts, control$options %||% list())
  res <- .Call(optcompR_ipopt_solve, as.numeric(x0), big(sf$lower),
               big(sf$upper), big(as.numeric(sf$expr_lower)),
               big(as.numeric(sf$expr_upper)), f, grad, g, jac,
               if (exact_h) hess else NULL, opts, environment(),
               as.integer(jidx[, 1L] - 1L), as.integer(jidx[, 2L] - 1L),
               as.integer(hidx[, 1L] - 1L), as.integer(hidx[, 2L] - 1L))
  term <- ipopt_termination(res$status)
  ok <- term %in% c("optimal", "almost_optimal")
  if (isTRUE(res$callback_error)) {
    why <- trimws(sub("^Error[^:]*: ", "", res$callback_message %||% ""))
    if (!ok) {
      stop(sprintf("An objective or constraint evaluation failed inside IPOPT: %s",
                   why), call. = FALSE)
    }
    warning(sprintf(paste("An evaluation failed at a trial point (%s);",
                          "IPOPT backtracked and converged."), why),
            call. = FALSE)
  }
  # at a limit, the returned point is a solution only if it is feasible
  tol <- 1e-6
  feas <- ok || (all(is.finite(res$x)) &&
    all(res$x >= sf$lower - tol * pmax(1, abs(sf$lower))) &&
    all(res$x <= sf$upper + tol * pmax(1, abs(sf$upper))) &&
    tryCatch({
      gx <- if (length(sf$expr_lower)) g(res$x) else numeric()
      all(gx >= sf$expr_lower - tol * pmax(1, abs(sf$expr_lower))) &&
        all(gx <= sf$expr_upper + tol * pmax(1, abs(sf$expr_upper)))
    }, error = function(e) FALSE))
  # IPOPT: grad(s f) + J' lambda - zL + zU = 0; optcompR: grad f = J'y + z
  list(termination = term,
       primal_status = if (feas) "feasible_point" else "infeasible_point",
       dual_status = if (ok) "feasible_point" else "not_available",
       x = res$x, objective = s * res$objective,
       row_dual = if (ok) -s * res$mult_g else NULL,
       col_dual = if (ok) s * (res$mult_x_L - res$mult_x_U) else NULL,
       iterations = NULL, gap = NA_real_, bound = NA_real_,
       message = sprintf("IPOPT status %d%s", res$status,
                         if (exact_h) "" else ", limited-memory Hessian"),
       local = cls$convexity != "convex", warm_start = any(!miss),
       extra = list(exact_hessian = exact_h), raw = res)
}

register_ipopt <- function() {
  if (!opt_ipopt_available()) return(invisible())
  register_solver(
    "ipopt", NULL,
    "IPOPT: interior point for nonlinear programs with exact Hessians (local)",
    classes = c("NLP", "QCQP", "QP", "LP"),
    global = character(),
    features = list(duals = TRUE, warm_start = TRUE, time_limit = TRUE,
                    iteration_limit = TRUE),
    solve = ipopt_solve_sf, version = function() "linked", priority = 65L,
    feasibility_tol = 1e-4)
}
