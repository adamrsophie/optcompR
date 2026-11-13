test_that("HiGHS uses start values and they invalidate the cache", {
  m <- opt_model("max")
  x <- opt_variable(m, "x", dim = 3, type = "integer", lower = 0, upper = 5)
  opt_constraint(m, sum(x) <= 7)
  opt_objective(m, sum(c(3, 2, 1) * x))
  k1 <- opt_standard_form(m)$key
  opt_set_start(m, x, c(5, 2, 0))
  expect_false(identical(opt_standard_form(m)$key, k1))
  expect_equal(opt_standard_form(m)$start, c(5, 2, 0))
  r <- opt_solve(m)
  expect_true(r$warm_start)
  expect_equal(opt_objective_value(r), 19)
  r2 <- opt_solve(m, warm_start = FALSE)
  expect_false(r2$warm_start)
})

test_that("creating models leaves the random number stream alone", {
  set.seed(11)
  a <- runif(1)
  set.seed(11)
  m <- opt_model()
  b <- runif(1)
  expect_identical(a, b)
})

test_that("user options override defaults without duplicates", {
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 1)
  opt_objective(m, x)
  # log_to_console is one of the defaults; changing threads instead would
  # ask HiGHS to resize its global scheduler, which it cannot do
  r <- opt_solve(m, options = list(log_to_console = FALSE, presolve = "off"))
  expect_identical(opt_status(r), "optimal")
  # a different thread count is refused before HiGHS sees it
  expect_error(opt_solve(m, options = list(threads = 3L)),
               "cannot change its thread count")
  expect_identical(opt_status(opt_solve(m, options = list(threads = 1L))),
                   "optimal")
})

test_that("MILP bounds include the objective constant", {
  m <- opt_model("max")
  x <- opt_variable(m, "x", type = "integer", lower = 0, upper = 4)
  opt_constraint(m, 2 * x <= 7)
  opt_objective(m, 5 * x + 100)
  r <- opt_solve(m)
  expect_equal(opt_objective_value(r), 115)
  expect_equal(opt_objective_bound(r), 115, tolerance = 1e-9)
  expect_equal(opt_objective_value(r, solver = TRUE), 115)
})

test_that("rows without variable terms and models without rows", {
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 0, upper = 2)
  p <- opt_parameter(m, "p", value = 1)
  opt_constraint(m, p <= 2) # becomes a row with no variables
  opt_objective(m, (x - 1.5)^2)
  r <- opt_solve(m)
  expect_equal(opt_value(r, x), 1.5, tolerance = 1e-6)
  skip_if_not_installed("scip")
  # SCIP works to absolute feasibility tolerance 1e-6 on the epigraph row,
  # so x is accurate to about its square root
  expect_equal(opt_value(opt_solve(m, solver = "scip"), x), 1.5,
               tolerance = 2e-3)
  skip_if_not_installed("osqp")
  m2 <- opt_model()
  y <- opt_variable(m2, "y")
  opt_objective(m2, (y - 3)^2)
  expect_equal(opt_value(opt_solve(m2, solver = "osqp"), y), 3,
               tolerance = 1e-5)
})

test_that("OSQP dual infeasibility is not called unbounded", {
  skip_if_not_installed("osqp")
  m <- opt_model()
  x <- opt_variable(m, "x")
  opt_objective(m, -x)
  expect_identical(opt_status(opt_solve(m, solver = "osqp")),
                   "infeasible_or_unbounded")
})

test_that("gradients at a result use the stored parameter values", {
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 1)
  p <- opt_parameter(m, "p", value = 3)
  opt_objective(m, p * x^2)
  r <- opt_solve(m)
  opt_set_parameter(m, p, 10)
  expect_equal(unname(opt_gradient(m, r)), 6, tolerance = 1e-8)
  expect_equal(unname(opt_gradient(m, c(x = 1))), 20)
})

test_that("double brackets do not expose internals", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 2)
  expect_error(x[[1]], "Use `\\[`")
})
