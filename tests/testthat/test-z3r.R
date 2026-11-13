test_that("Z3 decides feasibility exactly", {
  skip_if_not_installed("z3r")
  m <- opt_model()
  x <- opt_variable(m, "x", type = "integer", lower = 0, upper = 10)
  y <- opt_variable(m, "y", lower = 0)
  opt_constraint(m, 2 * x == 7)
  expect_identical(opt_z3_check(m)$status, "unsat")
  opt_delete_constraint(m, "c1")
  opt_constraint(m, 2 * x + 0.1 * y == 7)
  out <- opt_z3_check(m)
  expect_identical(out$status, "sat")
  expect_equal(2 * out$values[["x"]] + 0.1 * out$values[["y"]], 7)
  opt_objective(m, y)
  # y = 70 - 20x with integer x <= 3 gives y >= 10
  expect_identical(opt_z3_check(m, objective_bound = 10)$status, "sat")
  expect_identical(opt_z3_check(m, objective_bound = 9.99)$status, "unsat")
})

test_that("Z3 certifies solver optima and detects wrong claims", {
  skip_if_not_installed("z3r")
  m <- opt_model("max")
  x <- opt_variable(m, "x", type = "binary")
  y <- opt_variable(m, "y", type = "integer", lower = 0, upper = 3)
  opt_constraint(m, x + y <= 4)
  opt_objective(m, 3 * x + 2 * y)
  r <- opt_solve(m)
  expect_true(opt_z3_certify(r)$certified)
  # a falsified result: pretend the optimum were 8
  bad <- r
  bad$objective_value <- 8
  cert <- opt_z3_certify(bad)
  expect_false(cert$certified)
  expect_gt(3 * cert$better_point[["x"]] + 2 * cert$better_point[["y"]], 8)
})

test_that("nonlinear models are refused", {
  skip_if_not_installed("z3r")
  m <- opt_model()
  x <- opt_variable(m, "x")
  opt_constraint(m, x^2 <= 2)
  expect_error(opt_z3_check(m), "linear constraints only")
})
