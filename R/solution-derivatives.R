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
#' set. Mixed-integer and conic models are refused; for a mixed-integer
#' model, differentiate the [opt_fixed_relaxation()] instead.
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
  cls <- result$classification
  if (isTRUE(cls$integer)) {
    stop(paste("Solutions of mixed-integer models are not differentiable;",
               "differentiate opt_fixed_relaxation(result) instead."),
         call. = FALSE)
  }
  sf <- result$standard_form
  if (length(sf$cones)) {
    stop("Differentiation of conic solutions is not supported.", call. = FALSE)
  }
  model <- result$model
  ps <- resolve_symbols(model, parameter, "par")
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
  rs <- pmax(1, abs(ifelse(is.finite(lo), lo, 0)), abs(ifelse(is.finite(up), up, 0)))
  act <- (is.finite(lo) & abs(v - lo) <= sqrt(tol) * rs) |
    (is.finite(up) & abs(v - up) <= sqrt(tol) * rs)
  bs <- pmax(1, abs(ifelse(is.finite(sf$lower), sf$lower, 0)),
             abs(ifelse(is.finite(sf$upper), sf$upper, 0)))
  bact <- (is.finite(sf$lower) & abs(x - sf$lower) <= sqrt(tol) * bs) |
    (is.finite(sf$upper) & abs(x - sf$upper) <= sqrt(tol) * bs)
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
  J <- if (length(Jrows)) as.matrix(do.call(rbind, Jrows)) else matrix(0, 0, ns)
  H <- as.matrix(H)
  H <- (H + t(H)) / 2
  free <- cols[!bact]
  nf <- length(free)
  JA <- J[act, , drop = FALSE]
  na <- nrow(JA)
  K <- rbind(cbind(H[free, free, drop = FALSE], -t(JA[, free, drop = FALSE])),
             cbind(JA[, free, drop = FALSE], matrix(0, na, na)))
  rhs <- -rbind(H[free, ps, drop = FALSE], JA[, ps, drop = FALSE])
  sol <- tryCatch(solve(K, rhs), error = function(e) NULL)
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
  df <- gf[ps] - as.vector(crossprod(J[, ps, drop = FALSE], y))
  list(solution = dx, objective = stats::setNames(df, pn),
       active = c(sf$row_names[act],
                  if (any(bact)) paste0("bound:", sf$col_names[bact])),
       multipliers = multipliers)
}
