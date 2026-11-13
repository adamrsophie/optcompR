lp_model <- function(rhs = 4, cx = 3) {
  m <- opt_model("max")
  x <- opt_variable(m, "x", lower = 0)
  y <- opt_variable(m, "y", lower = 0)
  opt_constraint(m, x + y <= rhs, name = "cap")
  opt_constraint(m, x <= 2, name = "xmax")
  opt_constraint(m, y <= 3, name = "ymax")
  opt_objective(m, cx * x + 2 * y)
  m
}

test_that("ranging matches re-solves inside and at the range ends", {
  r <- opt_solve(lp_model())
  s <- opt_sensitivity(r)
  cap <- s$rows[s$rows$row == "cap", ]
  expect_equal(cap$dual, 2)
  # inside the rhs range the objective is linear with slope = dual
  mid <- (cap$bound_upper + 4) / 2
  expect_equal(opt_objective_value(opt_solve(lp_model(rhs = mid))),
               10 + 2 * (mid - 4), tolerance = 1e-8)
  expect_equal(opt_objective_value(opt_solve(lp_model(rhs = cap$bound_upper))),
               cap$objective_at_upper, tolerance = 1e-8)
  cx <- s$costs[s$costs$variable == "x", ]
  expect_equal(cx$lower, 2)
  # basis stays optimal: x stays 2 for coefficients inside the range
  r2 <- opt_solve(lp_model(cx = 2.5))
  expect_equal(opt_value(r2, opt_get(r2$model, "x")), 2, tolerance = 1e-8)
})

test_that("sensitivity is refused where it is not defined", {
  m <- opt_model()
  x <- opt_variable(m, "x", type = "integer", lower = 0, upper = 3)
  opt_objective(m, x)
  expect_error(opt_sensitivity(opt_solve(m)), "linear programs")
})

test_that("fixed relaxations give labeled conditional duals", {
  m <- opt_model("max")
  x <- opt_variable(m, "x", type = "integer", lower = 0, upper = 3)
  y <- opt_variable(m, "y", lower = 0)
  cap <- opt_constraint(m, x + y <= 4.5)
  opt_objective(m, 3 * x + 2 * y)
  r <- opt_solve(m)
  fr <- opt_fixed_relaxation(r)
  expect_true(fr$fixed_relaxation)
  expect_equal(unname(opt_dual(fr, cap)), 2)
  expect_identical(opt_variables(m)$type, c("integer", "continuous"))
  expect_identical(opt_classify(m)$class, "MILP")
})

test_that("scenario sweeps restore the parameter", {
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 0)
  d <- opt_parameter(m, "d", value = 1)
  opt_constraint(m, x >= d)
  opt_objective(m, 2 * x)
  sc <- opt_scenarios(m, d, c(1, 2, 3), keep = TRUE)
  expect_equal(sc$objective, c(2, 4, 6))
  expect_length(attr(sc, "results"), 3)
  expect_equal(m$value[optcompR:::resolve_symbols(m, d, "par")], 1)
})
