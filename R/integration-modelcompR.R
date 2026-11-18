# Least squares fitting of ODE models from modelcompR.

#' Least squares loss of an ODE model from modelcompR
#'
#' Builds the scalar expression
#' \eqn{\sum_{t, s} w_{ts} (u_s(t; \theta) - y_{ts})^2} for an
#' `ODEProblem` of the modelcompR package, with \eqn{\theta} given by
#' decision variables. The gradient is the derivative of the numerical
#' solution up to the integration tolerances. It comes from forward
#' sensitivities (`ODEForwardSensitivityProblem()`), which integrate states
#' times parameters equations, or from an adjoint through
#' `solve_gradient()` with `InterpolatingAdjoint()`, which integrates about
#' twice the states backward in time plus one quadrature per parameter. The
#' adjoint includes parameter-dependent initial values and events (such as
#' callbacks). Its derivative of the loss with respect to the saved states
#' is taken by central differences, which are exact for this quadratic loss
#' up to rounding. For a right-hand side given as an R function without a
#' Jacobian, modelcompR approximates the Jacobian products the
#' sensitivities need by finite differences. Value and gradient at the same
#' point share one solve.
#'
#' Observation times must lie in the integration interval. They may be in
#' any order and may repeat; data and weights rows follow their times, and
#' a state may be observed by several columns of `data`.
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
#'   parameters and more than twice as many parameters as states. A
#'   problem whose time origin is more than a million times its interval
#'   length always uses forward sensitivities. `states` index the full
#'   state vector; the problem's `save_idxs` is ignored.
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
  # one schedule for every method: the observation times must lie in the
  # integration interval, and they are sorted in the direction of
  # integration with data and weights; repeated times share a saved point
  span <- as.numeric(prob$tspan)
  lo <- min(span)
  hi <- max(span)
  slack <- 1e-12 * max(1, abs(span))
  if (any(!is.finite(times)) || any(times < lo - slack | times > hi + slack)) {
    stop(sprintf(paste("Observation times must lie in the integration",
                       "interval [%s, %s] of `prob`."), fmt_num(lo), fmt_num(hi)),
         call. = FALSE)
  }
  ut <- sort(unique(times), decreasing = span[2L] < span[1L])
  at <- match(times, ut)
  # the solves save the interval endpoints too, so the schedule is never a
  # single number (which modelcompR reads as a step size)
  grid <- unique(c(span[1L], ut, span[2L]))
  # `states` index the full state vector, so the internal solves save
  # every state whatever the problem's own save_idxs
  prob$kwargs$save_idxs <- NULL
  nu <- length(prob$u0)
  if (any(states < 1 | states > nu | states != round(states))) {
    stop(sprintf("`states` must be indices between 1 and %d.", nu), call. = FALSE)
  }
  np <- as_expr(theta)$len
  if (sensitivity == "auto") {
    # forward sensitivities integrate states x parameters equations, the
    # adjoint about twice the states plus one quadrature per parameter
    ns0 <- length(prob$u0)
    sensitivity <- if (np > 8L && np > 2L * ns0) "adjoint" else "forward"
  }
  # the backward solve loses the time step against a time origin far
  # larger than the interval (its gradient is wrong there), so such a
  # problem uses forward sensitivities
  if (sensitivity == "adjoint" &&
      max(abs(span)) > 1e6 * abs(span[2L] - span[1L])) {
    sensitivity <- "forward"
  }
  cache <- new.env(parent = emptyenv())
  saved_rows <- function(tt) {
    # the last saved value at each observation time (after an event
    # there): an exact time match first, else the nearest saved time
    # within rounding of it
    k <- vapply(ut, function(x) {
      k <- which(tt == x)
      if (!length(k)) {
        d <- abs(tt - x)
        k <- which(d == min(d) & d <= 4 * .Machine$double.eps * max(1, abs(x)))
      }
      if (!length(k)) NA_integer_ else k[length(k)]
    }, 1L)
    if (anyNA(k)) {
      stop("The ODE solve did not reach every observation time.",
           call. = FALSE)
    }
    k
  }
  states_at <- function(sol) {
    do.call(rbind, lapply(sol$u[saved_rows(as.numeric(sol$t))], as.numeric))
  }
  adjoint <- function(v) {
    p <- modelcompR::remake(prob, p = v)
    # the loss is differenced around the current states U0, written as
    # sum w ((U - U0)^2 + 2 (U0 - y)(U - U0)) plus the constant loss at U0,
    # so a large residual elsewhere does not swamp a small change
    U0 <- states_at(solve(p, alg, reltol = reltol, abstol = abstol,
                          saveat = grid))
    R0 <- U0[at, states, drop = FALSE] - data
    loss <- function(sol) {
      D <- states_at(sol)[at, states, drop = FALSE] -
        U0[at, states, drop = FALSE]
      sum(weights * (D^2 + 2 * R0 * D))
    }
    # modelcompR's solve_gradient composes the adjoint with the
    # initialization of parameter-dependent initial values and records
    # events for the backward solve
    g <- modelcompR::solve_gradient(loss, p, alg, saveat = grid,
                                    reltol = reltol, abstol = abstol,
                                    sensealg = modelcompR::InterpolatingAdjoint())
    list(value = sum(weights * R0^2), gradient = as.numeric(g$p))
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
                reltol = reltol, abstol = abstol, saveat = grid)
    s <- modelcompR::extract_local_sensitivities(fs, asmatrix = TRUE)
    k <- saved_rows(as.numeric(fs$t))[at]
    R <- s$u[k, states, drop = FALSE] - data
    g <- numeric(length(v))
    for (i in seq_along(times)) {
      S <- s$du[[k[i]]][states, , drop = FALSE]
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
