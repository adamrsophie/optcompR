fd_solution <- function(model, p, value, vars, h = 1e-5, ...) {
  at <- function(v) {
    opt_set_parameter(model, p, v)
    r <- opt_solve(model, ...)
    c(opt_value(r, vars), opt_objective_value(r))
  }
  out <- (at(value + h) - at(value - h)) / (2 * h)
  opt_set_parameter(model, p, value)
  out
}

test_that("QP solution derivatives match finite differences of re-solves", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 3, lower = c(-Inf, -Inf, 0))
  p <- opt_parameter(m, "p", value = 1)
  opt_constraint(m, sum(x) <= p, name = "budget")
  opt_constraint(m, x[1] - x[2] == 0.5 * p, name = "link")
  opt_objective(m, sum((x - c(2, 1, -1))^2) + p * x[1])
  r <- opt_solve(m)
  d <- opt_solution_jacobian(r, p)
  fd <- fd_solution(m, p, 1, x)
  expect_equal(as.vector(d$solution), fd[1:3], tolerance = 1e-5)
  expect_equal(unname(d$objective), fd[4], tolerance = 1e-5)
  expect_true("bound:x[3]" %in% d$active)
  expect_identical(d$multipliers, "solver")
})

test_that("LP derivatives with respect to right-hand sides and costs", {
  m <- opt_model("max")
  x <- opt_variable(m, "x", dim = 2, lower = 0)
  b <- opt_parameter(m, "b", value = c(4, 6))
  cst <- opt_parameter(m, "cst", value = 3)
  opt_constraint(m, x[1] + x[2] <= b[1])
  opt_constraint(m, x[1] + 3 * x[2] <= b[2])
  opt_objective(m, cst * x[1] + 4 * x[2])
  r <- opt_solve(m)
  d <- opt_solution_jacobian(r, b)
  expect_equal(dim(d$solution), c(2L, 2L))
  for (k in 1:2) {
    h <- 1e-4
    bv <- c(4, 6)
    opt_set_parameter(m, b, replace(bv, k, bv[k] + h))
    up <- opt_solve(m)
    opt_set_parameter(m, b, replace(bv, k, bv[k] - h))
    dn <- opt_solve(m)
    opt_set_parameter(m, b, bv)
    expect_equal(unname(d$solution[, k]),
                 (opt_value(up, x) - opt_value(dn, x)) / (2 * h), tolerance = 1e-6)
    expect_equal(unname(d$objective[k]),
                 (opt_objective_value(up) - opt_objective_value(dn)) / (2 * h),
                 tolerance = 1e-6)
  }
  # cost changes do not move a nondegenerate vertex
  dc <- opt_solution_jacobian(r, cst)
  expect_equal(as.vector(dc$solution), c(0, 0))
  expect_equal(unname(dc$objective), opt_value(r, x[1]))
})

test_that("nonlinear solution derivatives use exact second derivatives", {
  skip_if_not_installed("nloptr")
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 2, start = c(0.5, 0.5))
  p <- opt_parameter(m, "p", value = 0.3)
  opt_constraint(m, x[1]^2 + x[2]^2 <= 1 + p, name = "disk")
  opt_objective(m, exp(x[1] - p) - 2 * x[1] - 3 * x[2])
  r <- opt_solve(m)
  d <- opt_solution_jacobian(r, p)
  fd <- fd_solution(m, p, 0.3, x, h = 1e-4)
  expect_equal(as.vector(d$solution), fd[1:2], tolerance = 1e-4)
  expect_equal(unname(d$objective), fd[3], tolerance = 1e-4)
  expect_identical(d$multipliers, "estimated")
})

test_that("degenerate and discrete cases are refused or flagged", {
  m <- opt_model()
  x <- opt_variable(m, "x", type = "integer", lower = 0, upper = 5)
  p <- opt_parameter(m, "p", value = 2.5)
  opt_constraint(m, x >= p)
  opt_objective(m, x)
  expect_error(opt_solution_jacobian(opt_solve(m), p), "mixed-integer")
  # two identical active rows: LICQ fails
  m2 <- opt_model()
  y <- opt_variable(m2, "y", dim = 2)
  q <- opt_parameter(m2, "q", value = 1)
  opt_constraint(m2, y[1] + y[2] >= q)
  opt_constraint(m2, 2 * y[1] + 2 * y[2] >= 2 * q)
  opt_objective(m2, sum(y^2))
  expect_error(suppressWarnings(opt_solution_jacobian(opt_solve(m2), q)),
               "LICQ")
  # weakly active row: zero multiplier
  m3 <- opt_model()
  w <- opt_variable(m3, "w")
  s <- opt_parameter(m3, "s", value = 1)
  opt_constraint(m3, w <= s)
  opt_objective(m3, (w - s)^2)
  expect_warning(opt_solution_jacobian(opt_solve(m3), s), "Strict complementarity")
})
