test_that("IPOPT is compiled in where CI installs it", {
  skip_if_not(identical(Sys.getenv("OPTCOMPR_EXPECT_IPOPT"), "true"))
  expect_true(opt_ipopt_available())
})

skip_if_not(opt_ipopt_available(), "optcompR was built without IPOPT")

test_that("IPOPT solves nonlinear programs with exact derivatives", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 2, start = c(-1.2, 1))
  opt_objective(m, (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2)
  r <- opt_solve(m, solver = "ipopt")
  expect_identical(opt_status(r), "locally_optimal")
  expect_equal(opt_value(r, x), c(1, 1), tolerance = 1e-7)
  expect_true(r$extra$exact_hessian)
})

test_that("IPOPT duals are shadow prices", {
  build <- function(b, cap) {
    m <- opt_model("max")
    x <- opt_variable(m, "x", dim = 2, lower = 0)
    opt_constraint(m, x[1]^2 + x[2]^2 <= b, name = "disk")
    opt_constraint(m, x[1] <= cap, name = "cap")
    opt_objective(m, x[1] + x[2])
    m
  }
  r <- opt_solve(build(2, 0.5), solver = "ipopt")
  expect_identical(opt_status(r), "optimal")
  h <- 1e-5
  fd <- function(m1, m2) (opt_objective_value(opt_solve(m1, solver = "ipopt")) -
    opt_objective_value(opt_solve(m2, solver = "ipopt"))) / (2 * h)
  expect_equal(unname(opt_dual(r, "disk")),
               fd(build(2 + h, 0.5), build(2 - h, 0.5)), tolerance = 1e-5)
  expect_equal(unname(opt_dual(r, "cap")),
               fd(build(2, 0.5 + h), build(2, 0.5 - h)), tolerance = 1e-5)
  d <- opt_check(r)
  expect_true(d$feasible)
  expect_false(d$kkt$estimated)
  expect_true(d$kkt$pass)
})

test_that("IPOPT matches HiGHS on a convex QP, including bound duals", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 3, lower = c(0, 0, 1))
  opt_constraint(m, sum(x) == 2, name = "total")
  opt_objective(m, sum((x - c(1, -1, 0))^2))
  ri <- opt_solve(m, solver = "ipopt")
  rh <- opt_solve(m, solver = "highs")
  expect_equal(opt_value(ri, x), opt_value(rh, x), tolerance = 1e-7)
  expect_equal(unname(opt_dual(ri, "total")), unname(opt_dual(rh, "total")),
               tolerance = 1e-6)
  expect_equal(unname(ri$col_dual), unname(rh$col_dual), tolerance = 1e-6)
})

test_that("external functions without Hessians use limited memory", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 2, start = c(0, 0))
  f <- function(v) sum((v - c(1, 2))^2) + 0.1 * sum(v^4)
  g <- function(v) 2 * (v - c(1, 2)) + 0.4 * v^3
  opt_objective(m, opt_function(f, x, gradient = g))
  r <- opt_solve(m, solver = "ipopt")
  expect_false(r$extra$exact_hessian)
  expect_lt(max(abs(g(opt_value(r, x)))), 1e-6)
})

test_that("IPOPT reports infeasibility", {
  m <- opt_model()
  x <- opt_variable(m, "x")
  opt_constraint(m, x^2 <= -1)
  opt_objective(m, x)
  expect_identical(opt_status(opt_solve(m, solver = "ipopt")),
                   "locally_infeasible")
})
