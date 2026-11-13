# Reference problems. Independent references are built by hand and solved
# with highs::highs_solve() directly, never through optcompR's standard form.

direct <- function(...) highs::highs_solve(...)

test_that("LP reference: max 3x + 2y", {
  m <- opt_model("max")
  x <- opt_variable(m, "x", lower = 0)
  y <- opt_variable(m, "y", lower = 0)
  opt_constraint(m, x + y <= 4)
  opt_constraint(m, x <= 2)
  opt_constraint(m, y <= 3)
  opt_objective(m, 3 * x + 2 * y)
  r <- opt_solve(m)
  expect_identical(opt_status(r), "optimal")
  expect_equal(opt_value(r, x), 2, tolerance = 1e-8)
  expect_equal(opt_value(r, y), 2, tolerance = 1e-8)
  expect_equal(opt_objective_value(r), 10, tolerance = 1e-8)
  ref <- direct(L = c(3, 2), lower = c(0, 0), upper = c(Inf, Inf),
                A = rbind(c(1, 1), c(1, 0), c(0, 1)), lhs = rep(-Inf, 3),
                rhs = c(4, 2, 3), maximum = TRUE)
  expect_equal(opt_objective_value(r), ref$objective_value, tolerance = 1e-8)
  d <- opt_check(r)
  expect_true(d$feasible)
  expect_true(d$kkt$pass)
})

test_that("QP reference: min (x-1)^2 + (y-2)^2 s.t. x + y >= 4", {
  m <- opt_model()
  x <- opt_variable(m, "x")
  y <- opt_variable(m, "y")
  opt_constraint(m, x + y >= 4)
  opt_objective(m, (x - 1)^2 + (y - 2)^2)
  r <- opt_solve(m)
  expect_identical(opt_status(r), "optimal")
  expect_equal(opt_value(r, x), 1.5, tolerance = 1e-6)
  expect_equal(opt_value(r, y), 2.5, tolerance = 1e-6)
  expect_equal(opt_objective_value(r), 0.5, tolerance = 1e-6)
  # hand-built: 0.5 x'Qx + L'x + 5, Q = 2I, L = (-2, -4)
  ref <- direct(Q = 2 * diag(2), L = c(-2, -4), lower = c(-Inf, -Inf),
                upper = c(Inf, Inf), A = matrix(c(1, 1), 1), lhs = 4,
                rhs = Inf, offset = 5)
  expect_equal(opt_objective_value(r), ref$objective_value, tolerance = 1e-6)
  expect_true(opt_check(r)$kkt$pass)
})

test_that("MILP reference: max 3x + 2y, x binary, y integer", {
  m <- opt_model("max")
  x <- opt_variable(m, "x", type = "binary")
  y <- opt_variable(m, "y", type = "integer", lower = 0, upper = 3)
  opt_constraint(m, x + y <= 4)
  opt_objective(m, 3 * x + 2 * y)
  r <- opt_solve(m)
  expect_identical(opt_status(r), "optimal")
  expect_equal(c(opt_value(r, x), opt_value(r, y)), c(1, 3))
  expect_equal(opt_objective_value(r), 9)
  expect_error(opt_dual(r, "c1"), "mixed-integer")
  expect_identical(opt_dual_status(r), "not_available")
  ref <- direct(L = c(3, 2), lower = c(0, 0), upper = c(1, 3),
                A = matrix(c(1, 1), 1), lhs = -Inf, rhs = 4,
                types = c(2L, 2L), maximum = TRUE)
  expect_equal(opt_objective_value(r), ref$objective_value)
})

test_that("standard form matches a hand-built standard form", {
  m <- opt_model("max")
  x <- opt_variable(m, "x", dim = 3, lower = 0, upper = c(1, 2, 3))
  opt_constraint(m, x[1] + 2 * x[2] - x[3] <= 5)
  opt_constraint(m, 3 * x[1] + x[3] >= 1 + x[2])
  opt_constraint(m, x[1] - x[2], lower = -1, upper = 1)
  opt_objective(m, 4 * x[1] - x[2] + 2 * x[3] + 7)
  sf <- opt_standard_form(m)
  A <- rbind(c(1, 2, -1), c(3, -1, 1), c(1, -1, 0))
  expect_equal(as.matrix(sf$A), A, ignore_attr = TRUE)
  expect_equal(sf$row_lower, c(-Inf, 1, -1))
  expect_equal(sf$row_upper, c(5, Inf, 1))
  expect_equal(sf$c, c(4, -1, 2))
  expect_equal(sf$offset, 7)
  r <- opt_solve(m)
  ref <- direct(L = c(4, -1, 2), lower = c(0, 0, 0), upper = c(1, 2, 3),
                A = A, lhs = c(-Inf, 1, -1), rhs = c(5, Inf, 1),
                maximum = TRUE, offset = 7)
  expect_equal(opt_objective_value(r), ref$objective_value, tolerance = 1e-8)
})

test_that("reversed inequalities, signs and constants are not lost", {
  # x >= 3 must not become x <= 3
  m <- opt_model()
  x <- opt_variable(m, "x", lower = -10, upper = 10)
  opt_constraint(m, 3 <= x)
  opt_objective(m, x)
  expect_equal(opt_value(opt_solve(m), x), 3)
  # objective constants are reported
  opt_objective(m, x + 100)
  r <- opt_solve(m)
  expect_equal(opt_objective_value(r), 103)
  expect_equal(opt_objective_value(r, solver = TRUE), 103)
  # min vs max of the same expression
  opt_objective(m, 2 * x, sense = "max")
  expect_equal(opt_objective_value(opt_solve(m)), 20)
  opt_objective(m, 2 * x, sense = "min")
  expect_equal(opt_objective_value(opt_solve(m)), 6)
  # constants on both sides of a relation
  m2 <- opt_model("max")
  y <- opt_variable(m2, "y")
  opt_constraint(m2, 2 * y + 3 <= y + 10)
  opt_objective(m2, y)
  expect_equal(opt_value(opt_solve(m2), y), 7)
})

test_that("infeasible and unbounded models are reported", {
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 0)
  opt_constraint(m, x <= -1)
  opt_objective(m, x)
  r <- opt_solve(m)
  expect_identical(opt_status(r), "infeasible")
  expect_false(opt_has_solution(r))
  expect_error(opt_check(r), "no primal point")
  expect_warning(opt_value(r, x), "no primal solution")

  m <- opt_model("max")
  x <- opt_variable(m, "x", lower = 0)
  opt_objective(m, x)
  opt_constraint(m, x >= 1)
  r <- opt_solve(m)
  expect_true(opt_status(r) %in% c("unbounded", "infeasible_or_unbounded"))
  expect_false(isTRUE(opt_status(r) == "optimal"))

  m <- opt_model()
  x <- opt_variable(m, "x", type = "integer", lower = 0, upper = 10)
  opt_constraint(m, 2 * x == 3)
  opt_objective(m, x)
  expect_identical(opt_status(opt_solve(m)), "infeasible")
})

test_that("duals are shadow prices in both senses (finite differences)", {
  solve_obj <- function(sense, rhs, dir) {
    m <- opt_model(sense)
    x <- opt_variable(m, "x", lower = 0)
    y <- opt_variable(m, "y", lower = 0)
    cn <- if (dir == "le") opt_constraint(m, x + 2 * y <= rhs) else
      opt_constraint(m, x + 2 * y >= rhs)
    opt_constraint(m, x <= 3)
    if (sense == "max") opt_objective(m, 2 * x + 3 * y) else
      opt_objective(m, if (dir == "ge") 2 * x + 3 * y else -2 * x - 3 * y)
    r <- opt_solve(m)
    list(obj = opt_objective_value(r), dual = unname(opt_dual(r, cn)))
  }
  cases <- list(c("max", "le"), c("min", "ge"), c("min", "le"))
  for (cs in cases) {
    a <- solve_obj(cs[1], 4, cs[2])
    b <- solve_obj(cs[1], 4 + 1e-3, cs[2])
    expect_equal(a$dual, (b$obj - a$obj) / 1e-3, tolerance = 1e-6,
                 label = paste(cs, collapse = "/"))
  }
})

test_that("reduced costs follow the same convention", {
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 1)
  y <- opt_variable(m, "y", lower = 0)
  opt_constraint(m, x + y >= 0.5)
  opt_objective(m, 2 * x + y)
  r <- opt_solve(m)
  # raising the lower bound of x by d raises the objective by 2 d
  expect_equal(unname(opt_reduced_cost(r, x)), 2, tolerance = 1e-8)
})

test_that("degenerate, redundant and multiple-optima models solve", {
  m <- opt_model("max")
  x <- opt_variable(m, "x", lower = 0)
  y <- opt_variable(m, "y", lower = 0)
  opt_constraint(m, x + y <= 1)
  opt_constraint(m, 2 * x + 2 * y <= 2) # redundant
  opt_constraint(m, x <= 1)             # degenerate at the vertex (1, 0)
  opt_objective(m, x + y)               # every point of the edge is optimal
  r <- opt_solve(m)
  expect_identical(opt_status(r), "optimal")
  expect_equal(opt_objective_value(r), 1)
  expect_true(opt_check(r)$feasible)
})

test_that("fixed variables and feasibility problems", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 2, lower = 0, upper = 5)
  opt_fix(m, x[1], 2)
  opt_constraint(m, x[1] + x[2] == 3)
  r <- opt_solve(m) # no objective: feasibility problem
  expect_identical(opt_status(r), "optimal")
  expect_equal(opt_value(r, x), c(2, 1))
  expect_equal(opt_objective_value(r), 0)
})

test_that("large sparse LP solves and checks", {
  set.seed(42)
  n <- 2000
  k <- 500
  A <- Matrix::rsparsematrix(k, n, density = 0.005, rand.x = function(n) runif(n))
  b <- rep(10, k)
  cvec <- runif(n)
  m <- opt_model("max")
  x <- opt_variable(m, "x", dim = n, lower = 0, upper = 1)
  opt_constraint(m, A %*% x <= b)
  opt_objective(m, sum(cvec * x))
  r <- opt_solve(m)
  ref <- direct(L = cvec, lower = rep(0, n), upper = rep(1, n), A = A,
                lhs = rep(-Inf, k), rhs = b, maximum = TRUE)
  expect_equal(opt_objective_value(r), ref$objective_value, tolerance = 1e-7)
  expect_true(opt_check(r)$feasible)
})

test_that("MILP time limits do not claim optimality", {
  set.seed(7)
  n <- 120
  k <- 40
  A <- matrix(sample(1:40, n * k, TRUE), k)
  m <- opt_model("max")
  x <- opt_variable(m, "x", dim = n, type = "binary")
  opt_constraint(m, A %*% x <= rowSums(A) / 3)
  opt_objective(m, sum(sample(1:100, n, TRUE) * x))
  r <- opt_solve(m, time_limit = 0.05)
  if (opt_status(r) != "optimal") {
    expect_identical(opt_status(r), "time_limit")
    if (opt_has_solution(r)) expect_true(opt_check(r)$feasible)
  } else {
    expect_equal(opt_gap(r), 0, tolerance = 1e-4)
  }
})

test_that("name collisions are separate from constraint names", {
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 0)
  opt_constraint(m, x >= 1, name = "x")
  opt_objective(m, x)
  r <- opt_solve(m)
  expect_equal(unname(opt_dual(r, "x")), 1)
  expect_equal(opt_value(r, "x"), 1)
})
