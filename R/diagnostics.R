#' Independent feasibility and optimality diagnostics
#'
#' Recomputes, from the original model expressions and not from the solver
#' or the standard form, the constraint values, bound and integrality
#' violations and the objective value at the returned point. For continuous
#' problems with duals it also checks the Karush Kuhn Tucker conditions under
#' the sign convention documented in [opt_dual()]:
#'
#' * primal feasibility: rows and bounds hold within `tol`;
#' * stationarity: \eqn{\nabla f(x) - J(x)^\top y - z = 0}, with exact
#'   gradients;
#' * dual feasibility: each dual has the sign its active bound allows;
#' * complementary slackness: \eqn{y_i} times the row slack is zero.
#'
#' Violations are scaled by `max(1, |bound|)` for rows and by
#' `max(1, |gradient|)` for stationarity. Feasibility passes at `tol`; the
#' KKT test passes at `sqrt(tol)`, because duals and active sets carry the
#' square root of solver tolerances, and rows within `sqrt(tol)` of a bound
#' count as active. When a solver returns no multipliers, they are estimated
#' by least squares on the active set and the report says so. KKT conditions certify global
#' optimality only for convex problems; for nonconvex problems they describe
#' a stationary point.
#'
#' @param result An `opt_result`.
#' @param tol Feasibility and optimality tolerance.
#' @return An `opt_diagnostics` list with `max_row_violation`,
#'   `max_bound_violation`, `max_integrality_violation`, `objective_error`,
#'   `kkt` (or `NULL`), `feasible`, `rows`, a data frame of row residuals,
#'   and `special`, the violations of cone, indicator and SOS constraints
#'   (cones pass at `sqrt(tol)`). KKT conditions are not checked for conic
#'   models, whose multipliers live in the dual cones.
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x")
#' y <- opt_variable(m, "y")
#' opt_constraint(m, x + y >= 4)
#' opt_objective(m, (x - 1)^2 + (y - 2)^2)
#' opt_check(opt_solve(m))
#' @export
opt_check <- function(result, tol = 1e-6) {
  check_result(result)
  if (!opt_has_solution(result)) {
    stop(sprintf("The result has no primal point to check (status '%s').",
                 result$termination), call. = FALSE)
  }
  sf <- result$standard_form
  z <- result$z
  x <- z[sf$var_sym]
  # rows from the original expressions
  vals <- lapply(result$cons, function(b) eval_expr(b$expr, z))
  v <- unlist(vals, use.names = FALSE) %||% numeric()
  lo <- unlist(lapply(result$cons, `[[`, "lower"), use.names = FALSE) %||%
    numeric()
  up <- unlist(lapply(result$cons, `[[`, "upper"), use.names = FALSE) %||%
    numeric()
  viol <- pmax(lo - v, v - up, 0)
  scale <- pmax(1, abs(ifelse(is.finite(lo), lo, 0)),
                abs(ifelse(is.finite(up), up, 0)))
  rel <- viol / scale
  bviol <- pmax(sf$lower - x, x - sf$upper, 0) / pmax(1, abs(x))
  int <- sf$type %in% c("integer", "binary")
  iviol <- if (any(int)) abs(x[int] - round(x[int])) else 0
  obj_err <- if (is.na(result$solver_objective)) NA_real_ else
    abs(result$objective_value - result$solver_objective) /
      max(1, abs(result$objective_value))
  rows <- data.frame(row = sf$row_names, value = v, lower = lo, upper = up,
                     violation = viol, stringsAsFactors = FALSE)
  special <- special_violations(sf, x, tol)
  kkt <- NULL
  # conic KKT conditions involve cone multipliers and are not checked here
  if (!result$classification$integer && !length(sf$cones)) {
    if (!is.null(result$row_dual) && !is.null(result$col_dual)) {
      kkt <- kkt_check(result, v, lo, up, tol)
      kkt$estimated <- FALSE
    } else {
      est <- estimate_multipliers(result, v, lo, up, tol)
      kkt <- kkt_check(est, v, lo, up, tol)
      kkt$estimated <- TRUE
    }
  }
  sviol <- max(c(special$violation, 0))
  feasible <- max(c(rel, 0)) <= tol && max(c(bviol, 0)) <= tol &&
    max(c(iviol, 0)) <= tol && sviol <= sqrt(tol)
  structure(list(
    max_row_violation = max(c(rel, 0)),
    max_bound_violation = max(c(bviol, 0)),
    max_integrality_violation = max(c(iviol, 0)),
    max_special_violation = sviol, special = special,
    objective_error = obj_err, kkt = kkt, feasible = feasible,
    tol = tol, rows = rows, status = result$termination,
    convexity = result$classification$convexity
  ), class = "opt_diagnostics")
}

# Least-squares multipliers on the active set for solvers that return none.
estimate_multipliers <- function(result, v, lo, up, tol) {
  sf <- result$standard_form
  z <- result$z
  ns <- length(z)
  cols <- sf$var_sym
  x <- z[cols]
  n <- length(cols)
  g <- if (is.null(result$objective)) numeric(n) else
    as.vector(jac_expr(result$objective, z, ns)$J[, cols, drop = FALSE])
  rs <- pmax(1, abs(ifelse(is.finite(lo), lo, 0)), abs(ifelse(is.finite(up), up, 0)))
  act <- (is.finite(lo) & abs(v - lo) <= sqrt(tol) * rs) |
    (is.finite(up) & abs(v - up) <= sqrt(tol) * rs)
  bs <- pmax(1, abs(ifelse(is.finite(sf$lower), sf$lower, 0)),
             abs(ifelse(is.finite(sf$upper), sf$upper, 0)))
  bact <- (is.finite(sf$lower) & abs(x - sf$lower) <= sqrt(tol) * bs) |
    (is.finite(sf$upper) & abs(x - sf$upper) <= sqrt(tol) * bs)
  y <- numeric(length(v))
  zc <- numeric(n)
  if (any(act) || any(bact)) {
    J <- if (length(result$cons)) as.matrix(do.call(rbind, lapply(result$cons,
      function(b) jac_expr(b$expr, z, ns)$J[, cols, drop = FALSE]))) else
      matrix(0, 0, n)
    M <- cbind(t(J[act, , drop = FALSE]), diag(n)[, bact, drop = FALSE])
    sol <- qr.coef(qr(M), g)
    sol[is.na(sol)] <- 0
    # multiplier signs: s * y >= 0 at a lower bound, <= 0 at an upper bound,
    # free when both bounds hold (equality). A rank deficient system can
    # give an unrestricted solution with wrong signs although valid
    # multipliers exist, so refit with the signs enforced.
    s <- if (result$sense == "max") -1 else 1
    r_lo <- (is.finite(lo) & abs(v - lo) <= sqrt(tol) * rs)[act]
    r_up <- (is.finite(up) & abs(v - up) <= sqrt(tol) * rs)[act]
    b_lo <- (is.finite(sf$lower) & abs(x - sf$lower) <= sqrt(tol) * bs)[bact]
    b_up <- (is.finite(sf$upper) & abs(x - sf$upper) <= sqrt(tol) * bs)[bact]
    sg <- s * c(ifelse(r_lo & r_up, 0, ifelse(r_lo, 1, -1)),
                ifelse(b_lo & b_up, 0, ifelse(b_lo, 1, -1)))
    if (any(sg * sol < -1e-10 * max(1, abs(sol)))) {
      lb <- ifelse(sg > 0, 0, -Inf)
      ub <- ifelse(sg < 0, 0, Inf)
      fit <- stats::optim(pmin(pmax(sol, lb), ub),
        fn = function(u) sum((M %*% u - g)^2),
        gr = function(u) as.vector(2 * crossprod(M, M %*% u - g)),
        method = "L-BFGS-B", lower = lb, upper = ub,
        control = list(factr = 10, pgtol = 0, maxit = 1000L))
      sol <- fit$par
    }
    y[act] <- sol[seq_len(sum(act))]
    zc[bact] <- sol[sum(act) + seq_len(sum(bact))]
  }
  result$row_dual <- y
  result$col_dual <- zc
  result
}

kkt_check <- function(result, v, lo, up, tol) {
  sf <- result$standard_form
  z <- result$z
  ns <- length(z)
  cols <- sf$var_sym
  x <- z[cols]
  g <- if (is.null(result$objective)) numeric(length(cols)) else
    as.vector(jac_expr(result$objective, z, ns)$J[, cols, drop = FALSE])
  y <- unname(result$row_dual)
  zc <- unname(result$col_dual)
  Jty <- numeric(length(cols))
  if (length(result$cons)) {
    J <- do.call(rbind, lapply(result$cons, function(b) {
      jac_expr(b$expr, z, ns)$J[, cols, drop = FALSE]
    }))
    Jty <- as.vector(Matrix::crossprod(J, y))
  }
  stat <- g - Jty - zc
  stat_rel <- max(c(abs(stat) / pmax(1, abs(g)), 0))
  s <- if (result$sense == "max") -1 else 1
  # row dual signs: at upper bound y*s <= 0, at lower y*s >= 0
  rs <- pmax(1, abs(ifelse(is.finite(lo), lo, 0)), abs(ifelse(is.finite(up), up, 0)))
  at_lo <- is.finite(lo) & abs(v - lo) <= sqrt(tol) * rs
  at_up <- is.finite(up) & abs(v - up) <= sqrt(tol) * rs
  ys <- s * y
  dviol <- ifelse(at_lo & at_up, 0,
                  ifelse(at_up, pmax(ys, 0),
                         ifelse(at_lo, pmax(-ys, 0), abs(ys))))
  slack <- pmin(abs(v - lo), abs(up - v))
  comp <- abs(y) * ifelse(is.finite(slack), slack, 0)
  # bound duals
  bs <- pmax(1, abs(ifelse(is.finite(sf$lower), sf$lower, 0)),
             abs(ifelse(is.finite(sf$upper), sf$upper, 0)))
  bl <- is.finite(sf$lower) & abs(x - sf$lower) <= sqrt(tol) * bs
  bu <- is.finite(sf$upper) & abs(x - sf$upper) <= sqrt(tol) * bs
  zs <- s * zc
  zviol <- ifelse(bl & bu, 0, ifelse(bu, pmax(zs, 0),
                                     ifelse(bl, pmax(-zs, 0), abs(zs))))
  bslack <- pmin(abs(x - sf$lower), abs(sf$upper - x))
  bcomp <- abs(zc) * ifelse(is.finite(bslack), bslack, 0)
  list(stationarity = stat_rel,
       dual_feasibility = max(c(dviol, zviol, 0)) /
         max(1, abs(c(y, zc, 0))),
       complementarity = max(c(comp, bcomp, 0)),
       stationarity_residual = stat,
       pass = stat_rel <= sqrt(tol) &&
         max(c(dviol, zviol, 0)) <= sqrt(tol) * max(1, abs(c(y, zc, 0))) &&
         max(c(comp, bcomp, 0)) <= sqrt(tol))
}

#' @export
print.opt_diagnostics <- function(x, ...) {
  cat(sprintf("<opt_diagnostics: %s>\n", if (x$feasible) "feasible" else
    "INFEASIBLE point"))
  cat(sprintf("  max row violation (scaled):   %s\n", fmt_num(x$max_row_violation)))
  cat(sprintf("  max bound violation (scaled): %s\n", fmt_num(x$max_bound_violation)))
  cat(sprintf("  max integrality violation:    %s\n",
              fmt_num(x$max_integrality_violation)))
  if (nrow(x$special)) {
    cat(sprintf("  max cone/indicator/SOS violation: %s\n",
                fmt_num(x$max_special_violation)))
  }
  cat(sprintf("  objective recomputation error: %s\n", fmt_num(x$objective_error)))
  if (!is.null(x$kkt)) {
    cat(sprintf("  KKT: stationarity %s, dual feasibility %s, complementarity %s (%s)\n",
                fmt_num(x$kkt$stationarity), fmt_num(x$kkt$dual_feasibility),
                fmt_num(x$kkt$complementarity),
                if (x$kkt$pass) "pass" else "FAIL"))
    if (isTRUE(x$kkt$estimated)) {
      cat("  (multipliers estimated by least squares on the active set;",
          "the solver returned none)\n")
    }
    if (x$kkt$pass && x$convexity != "convex") {
      cat("  KKT conditions hold, but the problem is not known to be convex:",
          "this is a stationary point, not a certificate of global optimality.\n")
    }
  }
  invisible(x)
}
