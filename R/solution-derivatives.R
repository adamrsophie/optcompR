# Derivatives of optimal solutions with respect to parameters, by implicit
# differentiation of the KKT conditions on the active set.

#' Derivatives of the optimal solution with respect to parameters
#'
#' Differentiates the optimal point \eqn{x^*(p)} and optimal value
#' \eqn{f^*(p)} of a continuous model with respect to parameters, by the
#' implicit function theorem applied to the KKT conditions on the active
#' set at the solution. With \eqn{H} the Hessian of the Lagrangian
#' \eqn{L = f - y^\top g} in the free variables (those not at a bound),
#' \eqn{J_A} the Jacobian of the active rows and \eqn{y_A} their
#' multipliers, the derivatives solve
#' \deqn{\begin{pmatrix} H & -J_A^\top \\ J_A & 0 \end{pmatrix}
#'   \begin{pmatrix} dx/dp \\ dy_A/dp \end{pmatrix} =
#'   -\begin{pmatrix} \nabla^2_{xp} L \\ \partial g_A / \partial p
#'   \end{pmatrix}}
#' and \eqn{df^*/dp = \partial L / \partial p} (envelope theorem). All
#' second derivatives are exact.
#'
#' The result is valid when the solution is a strict local optimum that
#' satisfies linear independence of the active constraint gradients (LICQ),
#' the second-order sufficient condition and strict complementarity. The
#' function stops when the KKT matrix is singular (LICQ or second-order
#' sufficiency fails) and warns when an active row or bound has a zero
#' multiplier (strict complementarity fails, so the solution map may only
#' have one-sided derivatives). The active set is held fixed, so the
#' derivatives are local: a large parameter change can change the active
#' set. Mixed-integer models are refused; differentiate the
#' [opt_fixed_relaxation()] instead.
#'
#' Cone constraints are differentiated through their boundaries, which are
#' smooth away from the apex: an active second-order cone \eqn{\|u\| \le t}
#' becomes the row \eqn{\|u\| - t \le 0} (needs \eqn{u \ne 0}), an active
#' exponential cone the row \eqn{y e^{x/y} - z \le 0} (needs \eqn{y > 0}) and
#' an active power cone the row \eqn{|z| - x^lpha y^{1-lpha} \le 0}
#' (needs \eqn{z \ne 0} and \eqn{x, y > 0}); their multipliers are estimated
#' on the active set. An active PSD constraint \eqn{S(x, p) \succeq 0}
#' keeps the solver's multiplier \eqn{Z} (Clarabel and SCS return it) and
#' contributes the linearized complementarity
#' \eqn{\mathrm{sym}(dS\,Z + S\,dZ) = 0}, which is nonsingular under
#' strict complementarity and nondegeneracy. Inactive cones, including PSD
#' constraints whose matrix is positive definite, drop out. A cone at its
#' apex is refused.
#'
#' @param result A solved `opt_result` with a primal solution.
#' @param parameter A parameter expression or name.
#' @param tol Tolerance for active rows, bounds and zero multipliers.
#' @return A list with `solution` (matrix: variables by parameter elements),
#'   `objective` (gradient of the optimal value), `active` (names of active
#'   rows and bounds) and `multipliers` (`"solver"` or `"estimated"`).
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x", dim = 2)
#' p <- opt_parameter(m, "p", value = 1)
#' opt_constraint(m, sum(x) <= p, name = "budget")
#' opt_objective(m, sum((x - c(2, 1))^2))
#' r <- opt_solve(m)
#' opt_solution_jacobian(r, p)  # dx/dp = (0.5, 0.5)
#' @export
opt_solution_jacobian <- function(result, parameter, tol = 1e-7) {
  check_result(result)
  if (!opt_has_solution(result)) {
    stop("The result has no primal solution to differentiate.", call. = FALSE)
  }
  if (!result$termination %in% c("optimal", "locally_optimal",
                                 "almost_optimal")) {
    stop(sprintf(paste("The result terminated with `%s`, not at a feasible",
                       "optimum, so its point is not differentiable."),
                 result$termination), call. = FALSE)
  }
  cls <- result$classification
  if (isTRUE(cls$integer)) {
    stop(paste("Solutions of mixed-integer models are not differentiable;",
               "differentiate opt_fixed_relaxation(result) instead."),
         call. = FALSE)
  }
  if (length(result$standard_form$cones)) {
    result <- conic_as_smooth(result, tol)
  }
  sf <- result$standard_form
  model <- result$model
  ps <- result_symbols(result, parameter, "par")
  z <- result$z
  ns <- length(z)
  cols <- sf$var_sym
  x <- z[cols]
  n <- length(cols)
  multipliers <- "solver"
  # row values and bounds from the original expressions
  v <- unlist(lapply(result$cons, function(b) eval_expr(b$expr, z)),
              use.names = FALSE) %||% numeric()
  lo <- unlist(lapply(result$cons, `[[`, "lower"), use.names = FALSE) %||%
    numeric()
  up <- unlist(lapply(result$cons, `[[`, "upper"), use.names = FALSE) %||%
    numeric()
  if (is.null(result$row_dual) || is.null(result$col_dual)) {
    est <- estimate_multipliers(result, v, lo, up, tol)
    y <- est$row_dual
    zc <- est$col_dual
    multipliers <- "estimated"
  } else {
    y <- unname(result$row_dual)
    zc <- unname(result$col_dual)
  }
  chk <- result
  chk$row_dual <- y
  chk$col_dual <- zc
  kkt <- kkt_check(chk, v, lo, up, tol)
  if (!isTRUE(kkt$pass)) {
    stop(sprintf(paste("The point does not satisfy the KKT conditions",
                       "(stationarity %s, dual feasibility %s), so it is not a",
                       "differentiable optimum."),
                 fmt_num(kkt$stationarity), fmt_num(kkt$dual_feasibility)),
         call. = FALSE)
  }
  # primal feasibility, each side scaled by its own bound
  rviol <- bound_violation(v, lo, up)
  bviol <- bound_violation(x, sf$lower, sf$upper)
  if (max(c(rviol, bviol, 0)) > sqrt(tol)) {
    stop("The point violates the constraints, so it is not a differentiable optimum.",
         call. = FALSE)
  }
  nonsmooth_check(result, z, tol)
  # activity: tight within tol, or near (within sqrt(tol)) with a nonzero
  # multiplier; a near point with a zero multiplier is interior
  near <- function(d, b, m) {
    sc <- pmax(1, abs(b))
    is.finite(b) & (d <= tol * sc | (d <= sqrt(tol) * sc & abs(m) > sqrt(tol)))
  }
  act <- near(abs(v - lo), lo, y) | near(abs(v - up), up, y)
  bact <- near(abs(x - sf$lower), sf$lower, zc) |
    near(abs(x - sf$upper), sf$upper, zc)
  # strict complementarity concerns inequalities only
  ineq <- !(is.finite(lo) & lo == up)
  bineq <- !(sf$lower == sf$upper)
  weak <- c(sf$row_names[act & ineq & abs(y) <= tol],
            sf$col_names[bact & bineq & abs(zc) <= tol])
  if (length(weak)) {
    warning(sprintf(paste(
      "Strict complementarity fails at %s (active with a zero multiplier);",
      "the solution may only have one-sided derivatives."),
      paste(weak, collapse = ", ")), call. = FALSE)
  }
  # Lagrangian second derivatives over all symbols
  H <- if (is.null(result$objective)) sparse_zero(ns, ns) else
    hess_expr(result$objective, z, ns, 1)
  off <- 0L
  Jrows <- list()
  for (b in result$cons) {
    k <- b$expr$len
    w <- y[off + seq_len(k)]
    if (any(w != 0)) H <- H - hess_expr(b$expr, z, ns, w)
    Jrows <- c(Jrows, list(jac_expr(b$expr, z, ns)$J))
    off <- off + k
  }
  # H and J stay sparse over all symbols; only the blocks on the free
  # variables, active rows and parameters are made dense below
  J <- if (length(Jrows)) do.call(rbind, Jrows) else sparse_zero(0L, ns)
  psd <- result$extra$psd_active %||% list()
  s_sense <- if (result$sense == "max") -1 else 1
  for (q in psd) {
    H <- H - s_sense * hess_expr(q$vec, z, ns, svec_scale(q$k) * q$z)
  }
  H <- (H + Matrix::t(H)) / 2
  free <- cols[!bact]
  nf <- length(free)
  JA <- J[act, , drop = FALSE]
  na <- nrow(JA)
  Hff <- as.matrix(H[free, free, drop = FALSE])
  JAf <- as.matrix(JA[, free, drop = FALSE])
  # ponytail: the active-set KKT matrix is dense (free + active rows square);
  # a sparse LU with a condition estimate would be needed for huge active sets
  K <- rbind(cbind(Hff, -t(JAf)), cbind(JAf, matrix(0, na, na)))
  rhs <- -rbind(as.matrix(H[free, ps, drop = FALSE]),
                as.matrix(JA[, ps, drop = FALSE]))
  # active PSD constraints S(x, p) >= 0 with multiplier Z: the stationarity
  # term -s M' W dz and the linearized complementarity
  # sym(dS Z + S dZ) = 0, nonsingular under strict complementarity and
  # nondegeneracy (Alizadeh, Haeberly and Overton 1997)
  for (q in psd) {
    k <- q$k
    nz <- k * (k + 1L) / 2L
    w <- svec_scale(k)
    Mv <- as.matrix(jac_expr(q$vec, z, ns)$J)
    Zm <- smat_upper(q$z, k, scale = TRUE)
    rs <- qr(q$S, tol = sqrt(tol))$rank
    rz <- qr(Zm, tol = sqrt(tol))$rank
    if (rs + rz < k) {
      warning(sprintf(paste("Strict complementarity fails at the PSD",
                            "constraint '%s' (rank %d + %d < %d); the",
                            "solution may only have one-sided derivatives."),
                      q$name, rs, rz, k), call. = FALSE)
    }
    Cs <- -s_sense * t(Mv[, free, drop = FALSE]) %*% diag(w, nz)
    Cx <- vapply(free, function(j) svec_upper(smat_upper(Mv[, j], k) %*% Zm),
                 numeric(nz))
    Cz <- vapply(seq_len(nz), function(i) {
      e <- numeric(nz)
      e[i] <- 1
      svec_upper(q$S %*% smat_upper(e, k, scale = TRUE))
    }, numeric(nz))
    Rp <- -vapply(ps, function(j) svec_upper(smat_upper(Mv[, j], k) %*% Zm),
                  numeric(nz))
    nk <- nrow(K)
    pad <- matrix(0, nk - nf, nz)
    K <- rbind(cbind(K, rbind(Cs, pad)),
               cbind(matrix(Cx, nz, nf), matrix(0, nz, nk - nf),
                     matrix(Cz, nz, nz)))
    rhs <- rbind(rhs, matrix(Rp, nz, length(ps)))
  }
  bad <- ps[colSums(!is.finite(rhs)) > 0 |
              !is.finite(Matrix::colSums(J[, ps, drop = FALSE]))]
  if (length(bad)) {
    stop(sprintf(paste("The model is not differentiable in %s at this",
                       "solution (for example a parameter exponent of a",
                       "negative base)."),
                 paste(model$sym_name[bad], collapse = ", ")), call. = FALSE)
  }
  # second-order sufficiency: the Lagrangian Hessian (negated for a
  # maximum) is positive definite on the null space of the active rows
  # (with active PSD constraints the curvature condition has an extra
  # term, so nonsingularity of the KKT matrix is the test there)
  if (nf && !length(psd)) {
    A <- JAf
    rk <- if (na) qr(t(A))$rank else 0L
    if (rk < nf) {
      Z <- if (rk) qr.Q(qr(t(A)), complete = TRUE)[, (rk + 1L):nf, drop = FALSE]
        else diag(nf)
      s2 <- if (result$sense == "max") -1 else 1
      R <- s2 * crossprod(Z, Hff %*% Z)
      ev <- eigen((R + t(R)) / 2, symmetric = TRUE, only.values = TRUE)$values
      if (min(ev) <= sqrt(tol) * max(1, abs(ev))) {
        stop(paste(
          "The second-order sufficient condition fails: the Hessian of the",
          "Lagrangian is not positive definite on the null space of the",
          "active constraints, so the solution is not a strict local optimum",
          "and is not differentiable here."), call. = FALSE)
      }
    }
  }
  # every variable at a bound and no active row: the solution is locally
  # constant in the parameter
  sol <- if (!nrow(K)) matrix(0, 0, length(ps)) else
    tryCatch(solve(K, rhs), error = function(e) NULL)
  if (is.null(sol) || (length(K) && rcond(K) < 1e-12)) {
    stop(paste(
      "The KKT matrix on the active set is singular: the active constraint",
      "gradients are linearly dependent (LICQ fails) or the second-order",
      "sufficient condition fails, so the solution is not differentiable",
      "here."), call. = FALSE)
  }
  pn <- model$sym_name[ps]
  dx <- matrix(0, n, length(ps), dimnames = list(model$sym_name[cols], pn))
  dx[match(free, cols), ] <- sol[seq_len(nf), , drop = FALSE]
  gf <- if (is.null(result$objective)) numeric(ns) else
    as.vector(jac_expr(result$objective, z, ns)$J)
  df <- gf[ps] - as.vector(Matrix::crossprod(J[, ps, drop = FALSE], y))
  for (q in psd) {
    Jv <- as.matrix(jac_expr(q$vec, z, ns)$J)
    df <- df - s_sense * as.vector(crossprod(Jv[, ps, drop = FALSE],
                                             svec_scale(q$k) * q$z))
  }
  list(solution = dx, objective = stats::setNames(df, pn),
       active = c(sf$row_names[act],
                  if (any(bact)) paste0("bound:", sf$col_names[bact])),
       multipliers = multipliers)
}

# Smooth implicit differentiation is invalid at a kink: refuse when abs()
# is evaluated at (or within sqrt(tol) of) zero anywhere in the model.
nonsmooth_check <- function(result, z, tol) {
  kink <- function(e) {
    if (!is_expr(e) || e$type != "nl") return(FALSE)
    if (identical(e$op, "abs") &&
        any(abs(eval_expr(e$args[[1L]], z)) <= sqrt(tol))) return(TRUE)
    any(vapply(e$args %||% list(), kink, NA))
  }
  exprs <- c(list(result$objective), lapply(result$cons, `[[`, "expr"))
  if (any(vapply(exprs, kink, NA))) {
    stop(paste("abs() is evaluated at zero at this solution, where the model",
               "is not differentiable, so solution derivatives are not",
               "defined."), call. = FALSE)
  }
  invisible()
}

# Replace the cones of a conic result by equivalent smooth rows g <= 0 near
# the solution, so the NLP machinery applies. Inactive cones are dropped.
# The cones come from the result's own snapshot, so later edits of the model
# do not affect it.
conic_as_smooth <- function(result, tol) {
  z <- result$z
  near <- function(val, scale) abs(val) <= sqrt(tol) * max(1, scale)
  blocks <- list()
  psd <- list()
  for (it in result$special$cone %||% list()) {
    a <- lapply(it$args, as_expr)
    g <- switch(it$type,
      soc = {
        t <- a[[1L]]
        u <- a[[2L]]
        uv <- eval_expr(u, z)
        tv <- eval_expr(t, z)
        if (tv - sqrt(sum(uv^2)) > sqrt(tol) * max(1, abs(tv))) NULL else {
          if (near(sqrt(sum(uv^2)), abs(tv))) {
            stop(sprintf(paste("The second-order cone '%s' is active at its",
                               "apex, where it is not differentiable."),
                         it$name), call. = FALSE)
          }
          expr_unary("sqrt", expr_sum(expr_mul(u, u))) - t
        }
      },
      exp = {
        v <- vapply(a, eval_expr, 0, z = z)
        if (v[2L] <= sqrt(tol)) {
          stop(sprintf(paste("The exponential cone '%s' is at y = 0, where",
                             "it is not differentiable."), it$name),
               call. = FALSE)
        }
        gv <- v[2L] * exp(v[1L] / v[2L]) - v[3L]
        if (-gv > sqrt(tol) * max(1, abs(v[3L]))) NULL else
          a[[2L]] * exp(a[[1L]] / a[[2L]]) - a[[3L]]
      },
      pow = {
        v <- vapply(a, eval_expr, 0, z = z)
        al <- it$alpha
        cap <- max(v[1L], 0)^al * max(v[2L], 0)^(1 - al)
        if (cap - abs(v[3L]) > sqrt(tol) * max(1, cap)) NULL else {
          if (v[1L] <= sqrt(tol) || v[2L] <= sqrt(tol) || near(v[3L], cap)) {
            stop(sprintf(paste("The power cone '%s' is active where x, y or",
                               "z is zero, where it is not differentiable."),
                         it$name), call. = FALSE)
          }
          sign(v[3L]) * a[[3L]] - a[[1L]]^al * a[[2L]]^(1 - al)
        }
      },
      psd = {
        k <- it$k
        vec <- eval_expr(a[[1L]], z)
        Mv <- matrix(0, k, k)
        Mv[upper.tri(Mv, diag = TRUE)] <- vec
        Mv[lower.tri(Mv)] <- t(Mv)[lower.tri(Mv)]
        ev <- eigen(Mv, symmetric = TRUE, only.values = TRUE)$values
        if (min(ev) <= sqrt(tol) * max(1, abs(ev))) {
          zc <- result$extra[["cone_dual"]][[it$name]]
          if (is.null(zc)) {
            stop(sprintf(paste("The PSD constraint '%s' is active; its",
                               "derivatives need the cone multipliers of",
                               "Clarabel or SCS."), it$name), call. = FALSE)
          }
          psd <- c(psd, list(list(name = it$name, k = k, vec = a[[1L]],
                                  z = zc, S = Mv)))
        }
        NULL
      })
    if (!is.null(g)) {
      g$model <- result$model
      blocks <- c(blocks, list(list(name = it$name, expr = g, lower = -Inf,
                                    upper = 0, sense = "<=")))
    }
  }
  result$cons <- c(result$cons, blocks)
  sf <- result$standard_form
  sf$cones <- list()
  sf$row_names <- c(sf$row_names, vapply(blocks, `[[`, "", "name"))
  result$standard_form <- sf
  result$row_dual <- NULL
  result$col_dual <- NULL
  result$extra[["cone_dual"]] <- NULL
  # active PSD constraints keep their solver multipliers Z: their term
  # s * M' W z in the stationarity condition is a fixed gradient offset when
  # the other multipliers are estimated
  if (length(psd)) {
    s <- if (result$sense == "max") -1 else 1
    ns <- length(z)
    off <- numeric(ns)
    for (q in psd) {
      Jv <- jac_expr(q$vec, z, ns)$J
      off <- off + s * as.vector(Matrix::crossprod(Jv, svec_scale(q$k) * q$z))
    }
    result$extra$grad_offset <- off[sf$var_sym]
    result$extra$psd_active <- psd
  }
  result
}

# Symmetric matrix from upper-triangle-by-columns values, and back.
smat_upper <- function(v, k, scale = FALSE) {
  if (scale) v <- v / svec_scale(k)
  M <- matrix(0, k, k)
  M[upper.tri(M, diag = TRUE)] <- v
  M[lower.tri(M)] <- t(M)[lower.tri(M)]
  M
}
svec_upper <- function(M) {
  M <- (M + t(M)) / 2
  M[upper.tri(M, diag = TRUE)]
}
