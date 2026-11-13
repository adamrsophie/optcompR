test_that("box robust constraints use the worst case coefficients", {
  m <- opt_model("max")
  x <- opt_variable(m, "x", dim = 2, lower = 0)
  opt_robust_constraint(m, x, a0 = c(1, 1), b = 10, delta = c(0.5, 0.2))
  opt_objective(m, 3 * x[1] + 2 * x[2])
  r <- opt_solve(m)
  # worst case 1.5 x1 + 1.2 x2 <= 10
  expect_equal(opt_objective_value(r), 20, tolerance = 1e-7)
  expect_equal(opt_value(r, x), c(10 / 1.5, 0), tolerance = 1e-7)
  # x may be negative: the box uses |x|
  m2 <- opt_model("max")
  y <- opt_variable(m2, "y", lower = -5, upper = 5)
  opt_robust_constraint(m2, -y, a0 = 1, b = 2, delta = 1)
  opt_objective(m2, -y)
  expect_equal(opt_objective_value(opt_solve(m2)), 1, tolerance = 1e-7)
})

test_that("budget sets interpolate between nominal and box", {
  solve_budget <- function(gamma) {
    m <- opt_model("max")
    x <- opt_variable(m, "x", dim = 3, lower = 0)
    opt_robust_constraint(m, x, a0 = c(1, 1, 1), b = 6, delta = c(1, 1, 1),
                          uncertainty = "budget", gamma = gamma)
    opt_objective(m, sum(x))
    opt_objective_value(opt_solve(m))
  }
  expect_equal(solve_budget(0), 6, tolerance = 1e-7)
  expect_equal(solve_budget(3), 3, tolerance = 1e-7)
  # gamma = 1: by symmetry x = (t, t, t), 3t + t <= 6
  expect_equal(solve_budget(1), 4.5, tolerance = 1e-7)
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 2)
  expect_error(opt_robust_constraint(m, x, a0 = 1, b = 1, delta = 1,
                                     uncertainty = "budget", gamma = 3),
               "gamma")
})

test_that("ellipsoidal robust and Gaussian chance constraints are tight", {
  skip_if_not_installed("clarabel")
  P <- diag(c(0.5, 0.2))
  m <- opt_model("max")
  x <- opt_variable(m, "x", dim = 2, lower = 0)
  opt_robust_constraint(m, x, a0 = c(1, 1), b = 10, uncertainty = "ellipsoid",
                        P = P)
  opt_objective(m, 3 * x[1] + 2 * x[2])
  v <- opt_value(opt_solve(m), x)
  expect_equal(sum(v) + sqrt(sum((t(P) %*% v)^2)), 10, tolerance = 1e-6)
  # brute force over the boundary direction gives the same optimum
  th <- seq(0, pi / 2, length.out = 20001)
  d <- cbind(cos(th), sin(th))
  s <- 10 / (rowSums(d) + sqrt((0.5 * d[, 1])^2 + (0.2 * d[, 2])^2))
  expect_equal(3 * v[1] + 2 * v[2], max(s * (3 * d[, 1] + 2 * d[, 2])),
               tolerance = 1e-6)

  Sigma <- matrix(c(0.04, 0.01, 0.01, 0.02), 2)
  m <- opt_model("max")
  x <- opt_variable(m, "x", dim = 2, lower = 0)
  opt_chance_constraint(m, x, b = 10, prob = 0.95, mean = c(1, 1),
                        Sigma = Sigma)
  opt_objective(m, sum(x))
  v <- opt_value(opt_solve(m), x)
  p <- stats::pnorm((10 - sum(v)) / sqrt(drop(t(v) %*% Sigma %*% v)))
  expect_equal(p, 0.95, tolerance = 1e-6)
  expect_error(opt_chance_constraint(m, x, b = 1, prob = 0.3, mean = 1,
                                     Sigma = diag(2)), "not convex")
})

test_that("scenario chance constraints hold for every sample", {
  set.seed(1)
  S <- matrix(stats::rnorm(200, 1, 0.1), ncol = 2)
  m <- opt_model("max")
  x <- opt_variable(m, "x", dim = 2, lower = 0)
  opt_chance_constraint(m, x, b = 10, samples = S, method = "scenario")
  opt_objective(m, sum(x))
  v <- opt_value(opt_solve(m), x)
  expect_lte(max(S %*% v), 10 + 1e-7)
  expect_equal(max(S %*% v), 10, tolerance = 1e-7)
})

test_that("CVaR and worst case expectations match brute force", {
  R <- matrix(c(0.10, 0.02, -0.05, 0.04,
                0.03, 0.01, 0.02, 0.01), 4)
  grid <- seq(0, 1, by = 1e-4)
  loss_at <- function(w1) -(R %*% c(w1, 1 - w1))
  cvar_emp <- function(L, alpha, p) {
    # Rockafellar Uryasev minimum is attained at a scenario loss
    min(vapply(L, function(t) t + sum(p * pmax(L - t, 0)) / (1 - alpha), 0))
  }
  p <- c(0.1, 0.2, 0.3, 0.4)
  for (alpha in c(0.5, 0.9)) {
    m <- opt_model()
    w <- opt_variable(m, "w", dim = 2, lower = 0)
    opt_constraint(m, sum(w) == 1)
    opt_objective(m, opt_cvar(m, -(R %*% w), alpha = alpha, prob = p))
    r <- opt_solve(m)
    ref <- min(vapply(grid, function(g) cvar_emp(loss_at(g), alpha, p), 0))
    expect_equal(opt_objective_value(r), ref, tolerance = 1e-6)
    expect_equal(cvar_emp(loss_at(opt_value(r, w)[1]), alpha, p),
                 opt_objective_value(r), tolerance = 1e-7)
  }
  probs <- rbind(rep(0.25, 4), c(0.7, 0.1, 0.1, 0.1), c(0.1, 0.1, 0.7, 0.1))
  m <- opt_model()
  w <- opt_variable(m, "w", dim = 2, lower = 0)
  opt_constraint(m, sum(w) == 1)
  opt_objective(m, opt_worst_expectation(m, -(R %*% w), probs))
  ref <- min(vapply(grid, function(g) max(probs %*% loss_at(g)), 0))
  # the grid only bounds the optimum from above
  v <- opt_objective_value(opt_solve(m))
  expect_lte(v, ref + 1e-9)
  expect_lt(ref - v, 1e-5)
  expect_equal(as.numeric(optcompR:::eval_expr(opt_expectation(c(1, 3)), 0)), 2)
  expect_error(opt_cvar(m, -(R %*% w), alpha = 1), "alpha")
  expect_error(opt_cvar(m, -(R %*% w), alpha = 0.5, prob = c(1, 0)),
               "Probabilities")
})
