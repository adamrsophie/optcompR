# Least squares fitting of ODE models from modelcompR.

#' Least squares loss of an ODE model from modelcompR
#'
#' Builds the scalar expression
#' \eqn{\sum_{t, s} w_{ts} (u_s(t; \theta) - y_{ts})^2} for an
#' `ODEProblem` of the modelcompR package, with \eqn{\theta} given by
#' decision variables. The gradient is the exact derivative of the numerical
#' solution up to the integration tolerances; no finite differences are
#' used. It comes from forward sensitivities
#' (`ODEForwardSensitivityProblem()`), whose cost grows with the number of
#' parameters times the number of states, or from the adjoint equation
#' solved backward in time (`adjoint_sensitivities()`), whose cost does not
#' grow with the number of parameters. Value and gradient at the same point
#' share one solve.
#'
#' The loss is generally not convex, so solves are local: start from a
#' sensible `start` on the variables.
#'
#' @param prob A modelcompR `ODEProblem`; its parameters are replaced by
#'   `theta`.
#' @param theta Decision variables (or an expression) for the ODE
#'   parameters, in the order the right-hand side uses them.
#' @param times Observation times.
#' @param data Matrix of observations, one row per time and one column per
#'   observed state (`NA` for missing values).
#' @param states Indices of the observed states.
#' @param weights Matrix of weights like `data` (default 1).
#' @param reltol,abstol Integration tolerances.
#' @param alg Integration algorithm (default `modelcompR::Tsit5()`).
#' @param sensitivity `"forward"`, `"adjoint"`, or `"auto"` (the default),
#'   which uses the adjoint method when there are more than eight
#'   parameters and more than twice as many parameters as states.
#' @return A scalar nonlinear `opt_expr`.
#' @examplesIf requireNamespace("modelcompR", quietly = TRUE) && requireNamespace("nloptr", quietly = TRUE)
#' f <- function(u, p, t) -p[1] * u
#' prob <- modelcompR::ODEProblem(f, 1, c(0, 2), 1)
#' times <- c(0.5, 1, 1.5, 2)
#' y <- matrix(exp(-0.8 * times), ncol = 1)
#' m <- opt_model()
#' k <- opt_variable(m, "k", lower = 0, upper = 5, start = 0.3)
#' opt_objective(m, opt_ode_loss(prob, k, times, y))
#' opt_value(opt_solve(m), k)  # 0.8
#' @export
opt_ode_loss <- function(prob, theta, times, data, states = NULL,
                         weights = NULL, reltol = 1e-10, abstol = 1e-12,
                         alg = NULL,
                         sensitivity = c("auto", "forward", "adjoint")) {
  sensitivity <- match.arg(sensitivity)
  if (!requireNamespace("modelcompR", quietly = TRUE)) {
    stop("opt_ode_loss() needs the modelcompR package.", call. = FALSE)
  }
  check_numeric(times, "times")
  data <- as.matrix(data)
  if (nrow(data) != length(times)) {
    stop("`data` must have one row per observation time.", call. = FALSE)
  }
  states <- states %||% seq_len(ncol(data))
  if (length(states) != ncol(data)) {
    stop("`states` must have one entry per column of `data`.", call. = FALSE)
  }
  weights <- if (is.null(weights)) matrix(1, nrow(data), ncol(data)) else
    as.matrix(weights)
  if (!all(dim(weights) == dim(data))) {
    stop("`weights` must have the same shape as `data`.", call. = FALSE)
  }
  weights[is.na(data)] <- 0
  data[is.na(data)] <- 0
  alg <- alg %||% modelcompR::Tsit5()
  np <- as_expr(theta)$len
  if (sensitivity == "auto") {
    # forward sensitivities integrate states x parameters equations, the
    # adjoint about twice the states plus one quadrature per parameter
    ns0 <- length(prob$u0)
    sensitivity <- if (np > 8L && np > 2L * ns0) "adjoint" else "forward"
  }
  cache <- new.env(parent = emptyenv())
  adjoint <- function(v) {
    p <- modelcompR::remake(prob, p = v)
    sol <- solve(p, alg, reltol = reltol, abstol = abstol)
    U <- do.call(rbind, lapply(times, function(t) as.numeric(sol(t))))
    if (nrow(U) != length(times) || anyNA(U)) {
      stop("The ODE solve did not reach every observation time.",
           call. = FALSE)
    }
    R <- U[, states, drop = FALSE] - data
    nu <- ncol(U)
    res <- modelcompR::adjoint_sensitivities(
      sol, alg, t = times,
      dgdu_discrete = function(u, p, t, i) {
        g <- numeric(nu)
        g[states] <- 2 * weights[i, ] * (u[states] - data[i, ])
        g
      }, abstol = abstol, reltol = reltol)
    list(value = sum(weights * R^2), gradient = as.numeric(res$dp))
  }
  evaluate <- function(v) {
    if (identical(cache$at, v)) return(cache$out)
    if (sensitivity == "adjoint") {
      out <- adjoint(v)
      cache$at <- v
      cache$out <- out
      return(out)
    }
    p <- modelcompR::remake(prob, p = v)
    fs <- solve(modelcompR::ODEForwardSensitivityProblem(p), alg,
                reltol = reltol, abstol = abstol, saveat = times)
    s <- modelcompR::extract_local_sensitivities(fs, asmatrix = TRUE)
    U <- s$u[, states, drop = FALSE]
    if (nrow(U) != length(times)) {
      stop("The ODE solve did not reach every observation time.",
           call. = FALSE)
    }
    R <- U - data
    g <- numeric(length(v))
    for (i in seq_along(times)) {
      S <- s$du[[i]][states, , drop = FALSE]
      g <- g + 2 * as.vector(crossprod(S, weights[i, ] * R[i, ]))
    }
    out <- list(value = sum(weights * R^2), gradient = g)
    cache$at <- v
    cache$out <- out
    out
  }
  opt_function(function(v) evaluate(v)$value, theta,
               gradient = function(v) evaluate(v)$gradient, name = "ode_loss")
}
