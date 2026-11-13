enumerate_min <- function(f, grid) {
  vals <- apply(grid, 1, f)
  list(value = min(vals), at = unlist(grid[which.min(vals), ]))
}

test_that("MINLP branch and bound matches enumeration", {
  skip_if_not_installed("nloptr")
  m <- opt_model()
  x <- opt_variable(m, "x", type = "integer", lower = -3, upper = 3)
  y <- opt_variable(m, "y", type = "integer", lower = 0, upper = 4)
  u <- opt_variable(m, "u", lower = -5, upper = 5, start = 0)
  opt_constraint(m, x + y >= 2.5)
  opt_constraint(m, u >= x - 1)
  opt_objective(m, exp(0.5 * x) + (y - 1.7)^2 + u^2)
  cls <- opt_classify(m)
  expect_identical(cls$class, "MINLP")
  r <- opt_solve(m)
  expect_identical(r$solver, "nlpbb")
  # for fixed integers u* = max(x - 1, 0)
  grid <- expand.grid(x = -3:3, y = 0:4)
  grid <- grid[grid$x + grid$y >= 2.5, ]
  ref <- enumerate_min(function(g) exp(0.5 * g[1]) + (g[2] - 1.7)^2 +
                         max(g[1] - 1, 0)^2, grid)
  expect_equal(opt_objective_value(r), ref$value, tolerance = 1e-6)
  expect_equal(c(opt_value(r, x), opt_value(r, y)), unname(ref$at))
  expect_true(opt_check(r)$feasible)
})

test_that("convexity of a nonlinear relaxation is unknown, so results are local", {
  skip_if_not_installed("nloptr")
  m <- opt_model("max")
  x <- opt_variable(m, "x", type = "integer", lower = 0, upper = 6)
  opt_objective(m, sin(x))
  r <- opt_solve(m)
  expect_identical(opt_status(r), "locally_optimal")
  expect_true(is.na(opt_objective_bound(r)))
  expect_equal(opt_value(r, x), round(opt_value(r, x)))
})

test_that("infeasible MINLP and node limits", {
  skip_if_not_installed("nloptr")
  m <- opt_model()
  x <- opt_variable(m, "x", type = "integer", lower = 0, upper = 3)
  opt_constraint(m, 2 * x == 3)
  opt_objective(m, exp(x))
  expect_true(opt_status(opt_solve(m, solver = "nlpbb")) %in%
                c("infeasible", "locally_infeasible"))
  m2 <- opt_model()
  z <- opt_variable(m2, "z", dim = 6, type = "integer", lower = 0, upper = 9)
  opt_objective(m2, sum(exp(0.3 * z)) + sum((z - 4.5)^2))
  r <- opt_solve(m2, options = list(max_nodes = 3L))
  expect_true(opt_status(r) %in% c("node_limit"))
})

test_that("indicator and SOS constraints hold on the nonlinear path", {
  skip_if_not_installed("nloptr")
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 0, upper = 8, start = 1)
  on <- opt_variable(m, "on", type = "binary")
  opt_indicator(m, on, x <= 0, active = 0)
  opt_constraint(m, x >= 3)
  opt_objective(m, exp(0.2 * x) + 10 * on)
  # without the indicator, x = 3 and on = 0 would be cheaper
  r <- opt_solve(m, solver = "nlpbb")
  expect_equal(opt_value(r, c(x, on)), c(3, 1), tolerance = 1e-6)
  expect_true(opt_check(r)$feasible)
  expect_equal(max(opt_check(r)$special$violation), 0)
  m2 <- opt_model()
  v <- opt_variable(m2, "v", dim = 3, lower = 0, upper = 2, start = c(1, 1, 1))
  opt_sos1(m2, v)
  opt_objective(m2, sum((v - c(1.5, 1, 0.5))^2) + exp(-sum(v)))
  r2 <- opt_solve(m2, solver = "nlpbb")
  vv <- opt_value(r2, v)
  expect_lte(sum(abs(vv) > 1e-6), 1)
  # enumerate the three single-support candidates
  f1 <- function(i) {
    o <- stats::optimize(function(t) {
      w <- replace(numeric(3), i, t)
      sum((w - c(1.5, 1, 0.5))^2) + exp(-sum(w))
    }, c(0, 2), tol = 1e-10)
    o$objective
  }
  expect_equal(opt_objective_value(r2), min(vapply(1:3, f1, 0)),
               tolerance = 1e-6)
})
