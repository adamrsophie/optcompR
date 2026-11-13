# Robust and stochastic optimization through exact reformulations.

aux_name <- function(model, base) {
  k <- 1L
  repeat {
    nm <- sprintf(".%s%d", base, k)
    if (is.na(index_get(model$block_index, nm))) return(nm)
    k <- k + 1L
  }
}

affine_vector <- function(model, x, what = "x") {
  x <- check_affine(model_expr(model, x, what), what)
  x
}

#' Robust linear constraints
#'
#' Imposes \eqn{a^\top x \le b} for every coefficient vector \eqn{a} in an
#' uncertainty set, through an exact reformulation:
#'
#' * `"box"`: \eqn{a_i \in [a_{0i} - \delta_i, a_{0i} + \delta_i]}, giving
#'   \eqn{a_0^\top x + \delta^\top |x| \le b} (linear, with auxiliary
#'   variables for \eqn{|x|}).
#' * `"ellipsoid"`: \eqn{a = a_0 + P u} with \eqn{\|u\|_2 \le 1}, giving
#'   \eqn{a_0^\top x + \|P^\top x\|_2 \le b} (a second-order cone).
#' * `"budget"`: box deviations of which at most `gamma` act at once
#'   (Bertsimas and Sim), reformulated by linear duality.
#'
#' @param model An [opt_model()].
#' @param x Affine vector expression, usually decision variables.
#' @param a0 Nominal coefficients.
#' @param b Right-hand side.
#' @param uncertainty `"box"`, `"ellipsoid"` or `"budget"`.
#' @param delta Half-widths of the box (box and budget sets).
#' @param P Matrix with `length(x)` rows shaping the ellipsoid.
#' @param gamma Budget of uncertainty, between 0 and `length(x)`.
#' @param name Base name of the generated constraints.
#' @return A list of the generated constraint references, invisibly.
#' @examples
#' m <- opt_model("max")
#' x <- opt_variable(m, "x", dim = 2, lower = 0)
#' opt_robust_constraint(m, x, a0 = c(1, 1), b = 10, delta = c(0.5, 0.2))
#' opt_objective(m, 3 * x[1] + 2 * x[2])
#' opt_value(opt_solve(m), x)
#' @export
opt_robust_constraint <- function(model, x, a0, b,
                                  uncertainty = c("box", "ellipsoid", "budget"),
                                  delta = NULL, P = NULL, gamma = NULL,
                                  name = NULL) {
  check_model(model)
  uncertainty <- match.arg(uncertainty)
  x <- affine_vector(model, x)
  n <- x$len
  check_numeric(a0, "a0")
  a0 <- recycle_to(a0, n, "a0")
  check_numeric(b, "b")
  if (length(b) != 1L) stop("`b` must be a number.", call. = FALSE)
  base <- name %||% aux_name(model, "robust")
  nominal <- sum(a0 * x)
  if (uncertainty == "ellipsoid") {
    if (is.null(P)) stop("The ellipsoid needs `P`.", call. = FALSE)
    P <- as.matrix(P)
    if (nrow(P) != n) stop("`P` must have one row per element of `x`.",
                           call. = FALSE)
    return(invisible(list(opt_soc_constraint(model, t(P) %*% x, b - nominal,
                                             name = base))))
  }
  if (is.null(delta)) stop("Box and budget sets need `delta`.", call. = FALSE)
  check_numeric(delta, "delta")
  delta <- recycle_to(delta, n, "delta")
  if (any(delta < 0)) stop("`delta` must be nonnegative.", call. = FALSE)
  s <- opt_variable(model, paste0(base, ".abs"), dim = n, lower = 0)
  refs <- list(opt_constraint(model, s - x >= 0, name = paste0(base, ".pos")),
               opt_constraint(model, s + x >= 0, name = paste0(base, ".neg")))
  if (uncertainty == "box") {
    refs <- c(refs, list(opt_constraint(model, nominal + sum(delta * s) <= b,
                                        name = base)))
    return(invisible(refs))
  }
  check_numeric(gamma, "gamma")
  if (length(gamma) != 1L || gamma < 0 || gamma > n) {
    stop("`gamma` must be between 0 and the length of `x`.", call. = FALSE)
  }
  # max over |S| <= gamma (fractional allowed) of sum delta_i s_i equals
  # min gamma p + sum q  s.t.  p + q_i >= delta_i s_i,  p, q >= 0
  p <- opt_variable(model, paste0(base, ".p"), lower = 0)
  q <- opt_variable(model, paste0(base, ".q"), dim = n, lower = 0)
  refs <- c(refs, list(
    opt_constraint(model, rep(p, n) + q - delta * s >= 0,
                   name = paste0(base, ".dual")),
    opt_constraint(model, nominal + gamma * p + sum(q) <= b, name = base)))
  invisible(refs)
}

#' Chance constraints
#'
#' Requires \eqn{P(a^\top x \le b) \ge p} for a random coefficient vector
#' \eqn{a}.
#'
#' * `method = "gaussian"`: \eqn{a \sim N(\mu, \Sigma)}. The constraint is
#'   equivalent to \eqn{\mu^\top x + \Phi^{-1}(p) \|\Sigma^{1/2} x\|_2 \le b},
#'   a second-order cone for \eqn{p \ge 0.5}. This is exact.
#' * `method = "scenario"`: the constraint is imposed for every row of
#'   `samples`. This is the scenario approximation: with \eqn{N} independent
#'   samples and \eqn{d} decision variables, the solution violates the chance
#'   constraint with probability at most \eqn{\epsilon} when
#'   \eqn{N \ge \frac{2}{\epsilon}(\ln\frac{1}{\beta} + d)} (Calafiore and
#'   Campi), with confidence \eqn{1 - \beta}. It is labeled an
#'   approximation, not an exact reformulation.
#'
#' @inheritParams opt_robust_constraint
#' @param prob Required probability.
#' @param mean,Sigma Mean vector and covariance matrix (Gaussian method).
#' @param samples Matrix of sampled coefficient vectors, one per row
#'   (scenario method).
#' @param method `"gaussian"` or `"scenario"`.
#' @return The constraint reference, invisibly.
#' @examples
#' m <- opt_model("max")
#' x <- opt_variable(m, "x", dim = 2, lower = 0)
#' opt_chance_constraint(m, x, b = 10, prob = 0.95, mean = c(1, 1),
#'                       Sigma = diag(c(0.04, 0.01)))
#' opt_objective(m, sum(x))
#' if (requireNamespace("clarabel", quietly = TRUE)) opt_value(opt_solve(m), x)
#' @export
opt_chance_constraint <- function(model, x, b, prob, mean = NULL, Sigma = NULL,
                                  samples = NULL,
                                  method = c("gaussian", "scenario"),
                                  name = NULL) {
  check_model(model)
  method <- match.arg(method)
  x <- affine_vector(model, x)
  n <- x$len
  if (method == "scenario") {
    if (is.null(samples)) stop("The scenario method needs `samples`.",
                               call. = FALSE)
    S <- as.matrix(samples)
    if (ncol(S) != n) stop("`samples` must have one column per element of `x`.",
                           call. = FALSE)
    ref <- opt_constraint(model, S %*% x <= b, name = name)
    return(invisible(ref))
  }
  if (!is.numeric(prob) || length(prob) != 1L || prob < 0.5 || prob >= 1) {
    stop(paste("`prob` must be in [0.5, 1): below 0.5 the Gaussian chance",
               "constraint is not convex."), call. = FALSE)
  }
  mean <- recycle_to(as.numeric(mean), n, "mean")
  Sigma <- as.matrix(Sigma)
  if (!all(dim(Sigma) == n) || !isSymmetric(Sigma)) {
    stop("`Sigma` must be a symmetric matrix of size length(x).", call. = FALSE)
  }
  ev <- eigen(Sigma, symmetric = TRUE)
  if (any(ev$values < -1e-10 * max(1, abs(ev$values)))) {
    stop("`Sigma` must be positive semidefinite.", call. = FALSE)
  }
  R <- diag(sqrt(pmax(ev$values, 0)), n) %*% t(ev$vectors)
  kq <- stats::qnorm(prob)
  rhs <- b - sum(mean * x)
  if (kq == 0) return(invisible(opt_constraint(model, rhs >= 0, name = name)))
  invisible(opt_soc_constraint(model, kq * (R %*% x), rhs, name = name))
}

#' Risk measures and expectations over scenarios
#'
#' * `opt_expectation(losses, prob)`: the sample average (or weighted
#'   expectation) of scenario losses, the objective of sample average
#'   approximation.
#' * `opt_cvar(model, losses, alpha, prob)`: conditional value at risk at
#'   level `alpha`, through the Rockafellar and Uryasev formulation
#'   \eqn{t + \frac{1}{1-\alpha} \sum_k p_k (L_k - t)_+} with auxiliary
#'   variables. Minimizing it, or bounding it from above, is exact.
#' * `opt_worst_expectation(model, losses, probs)`: distributionally robust
#'   expectation \eqn{\max_j p_j^\top L} over a finite ambiguity set of
#'   distributions (rows of `probs`), and so over its convex hull, by an
#'   epigraph variable. Minimizing it, or bounding it from above, is exact.
#'
#' The second and third return expressions built on new variables and
#' constraints; they model the risk measure only in the direction stated
#' (in a minimization objective or on the left of `<=`).
#'
#' @param model An [opt_model()].
#' @param losses Expression with one element per scenario.
#' @param alpha Confidence level in (0, 1).
#' @param prob Scenario probabilities (default uniform).
#' @param probs Matrix of candidate distributions, one per row.
#' @param name Base name for generated variables.
#' @return A scalar expression.
#' @examples
#' # portfolio over 4 return scenarios, minimizing CVaR of the loss
#' R <- matrix(c(0.10, 0.02, -0.05, 0.04,
#'               0.03, 0.01,  0.02, 0.01), 4)
#' m <- opt_model()
#' w <- opt_variable(m, "w", dim = 2, lower = 0)
#' opt_constraint(m, sum(w) == 1)
#' loss <- -(R %*% w)
#' opt_objective(m, opt_cvar(m, loss, alpha = 0.75))
#' opt_value(opt_solve(m), w)
#' @export
opt_cvar <- function(model, losses, alpha, prob = NULL, name = NULL) {
  check_model(model)
  L <- model_expr(model, losses, "losses")
  K <- L$len
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be in (0, 1).", call. = FALSE)
  }
  prob <- scenario_prob(prob, K)
  base <- name %||% aux_name(model, "cvar")
  t <- opt_variable(model, paste0(base, ".t"))
  u <- opt_variable(model, paste0(base, ".excess"), dim = K, lower = 0)
  opt_constraint(model, u - L + rep(t, K) >= 0, name = paste0(base, ".excess"))
  t + sum(prob * u) / (1 - alpha)
}

#' @rdname opt_cvar
#' @export
opt_expectation <- function(losses, prob = NULL) {
  L <- as_expr(losses)
  sum(scenario_prob(prob, L$len) * L)
}

#' @rdname opt_cvar
#' @export
opt_worst_expectation <- function(model, losses, probs, name = NULL) {
  check_model(model)
  L <- model_expr(model, losses, "losses")
  probs <- as.matrix(probs)
  if (ncol(probs) != L$len) stop("`probs` must have one column per scenario.",
                                 call. = FALSE)
  for (j in seq_len(nrow(probs))) scenario_prob(probs[j, ], L$len)
  base <- name %||% aux_name(model, "worst")
  w <- opt_variable(model, paste0(base, ".w"))
  opt_constraint(model, rep(w, nrow(probs)) - probs %*% L >= 0,
                 name = paste0(base, ".epigraph"))
  w
}

scenario_prob <- function(prob, K) {
  if (is.null(prob)) return(rep(1 / K, K))
  check_numeric(prob, "prob")
  if (length(prob) != K || any(prob < 0) || abs(sum(prob) - 1) > 1e-9) {
    stop("Probabilities must be nonnegative, one per scenario, summing to 1.",
         call. = FALSE)
  }
  prob
}
