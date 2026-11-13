# Exact derivatives of model functions.

point_vector <- function(model, at) {
  vs <- var_symbols(model)
  if (inherits(at, "opt_result")) {
    if (!identical(at$model, model)) stop("`at` is a result of another model.",
                                          call. = FALSE)
    at <- at$z[vs]
  }
  check_numeric(at, "at")
  if (!is.null(names(at)) && all(nzchar(names(at)))) {
    nm <- model$sym_name[vs]
    if (!setequal(names(at), nm)) {
      stop("Names of `at` must match the variable names of the model.",
           call. = FALSE)
    }
    at <- at[nm]
  }
  if (length(at) != length(vs)) {
    stop(sprintf("`at` must have one value per variable (%d).", length(vs)),
         call. = FALSE)
  }
  z <- model$value
  z[is.na(z)] <- 0
  z[vs] <- as.numeric(at)
  list(z = z, vs = vs)
}

#' Exact derivatives of the objective and constraints
#'
#' `opt_gradient()` returns the gradient of the objective,
#' `opt_jacobian()` the sparse Jacobian of all constraint rows (row order of
#' [opt_standard_form()]), and `opt_hessian()` the Hessian of the objective,
#' or, when `multipliers` are given, the Hessian of the Lagrangian
#' \eqn{L = f(x) - \sum_i y_i g_i(x)} in the sign convention of
#' [opt_dual()]. Derivatives are exact: polynomial parts are differentiated
#' structurally (for \eqn{f = \frac12 x^\top Q x + c^\top x} the gradient is
#' \eqn{\frac12 (Q + Q^\top) x + c}), and nonlinear nodes use the chain rule
#' with sparse forward-mode Jacobians and exact second derivatives. No finite
#' differences are used. Parameters take their current values. Functions
#' that are not differentiable at a point (for example `abs()` at zero)
#' return a one-sided value; check the point when it matters.
#'
#' @param model An [opt_model()].
#' @param at Point: a numeric vector with one value per variable (named by
#'   variable names, or in model order), or an `opt_result`.
#' @param multipliers Optional row multipliers for the Lagrangian Hessian.
#' @return A named numeric vector (gradient) or a sparse `Matrix` (Jacobian,
#'   Hessian) with dimnames.
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x")
#' y <- opt_variable(m, "y")
#' opt_objective(m, x^2 + 3 * x * y + exp(y))
#' opt_constraint(m, x * y <= 4)
#' opt_gradient(m, c(1, 2))
#' opt_jacobian(m, c(1, 2))
#' opt_hessian(m, c(1, 2))
#' @export
opt_gradient <- function(model, at) {
  check_model(model)
  p <- point_vector(model, at)
  ns <- length(p$z)
  if (is.null(model$objective)) {
    g <- numeric(length(p$vs))
  } else {
    g <- as.vector(jac_expr(model$objective, p$z, ns)$J[, p$vs, drop = FALSE])
  }
  stats::setNames(g, model$sym_name[p$vs])
}

#' @rdname opt_gradient
#' @export
opt_jacobian <- function(model, at) {
  check_model(model)
  p <- point_vector(model, at)
  ns <- length(p$z)
  if (length(model$cons) == 0L) {
    return(sparse_zero(0L, length(p$vs)))
  }
  J <- do.call(rbind, lapply(model$cons, function(b) {
    jac_expr(b$expr, p$z, ns)$J[, p$vs, drop = FALSE]
  }))
  dimnames(J) <- list(opt_standard_form(model)$row_names, model$sym_name[p$vs])
  J
}

#' @rdname opt_gradient
#' @export
opt_hessian <- function(model, at, multipliers = NULL) {
  check_model(model)
  p <- point_vector(model, at)
  ns <- length(p$z)
  H <- sparse_zero(ns, ns)
  if (!is.null(model$objective)) H <- hess_expr(model$objective, p$z, ns, 1)
  if (!is.null(multipliers)) {
    nr <- sum(vapply(model$cons, function(b) b$expr$len, 0L))
    multipliers <- recycle_to(as.numeric(multipliers), nr, "multipliers")
    off <- 0L
    for (b in model$cons) {
      k <- b$expr$len
      w <- multipliers[off + seq_len(k)]
      if (any(w != 0)) H <- H - hess_expr(b$expr, p$z, ns, w)
      off <- off + k
    }
  }
  H <- H[p$vs, p$vs, drop = FALSE]
  H <- Matrix::forceSymmetric(methods::as(H, "generalMatrix"))
  dimnames(H) <- list(model$sym_name[p$vs], model$sym_name[p$vs])
  H
}
