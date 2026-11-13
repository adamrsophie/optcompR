rosen <- function(v) (1 - v[1])^2 + 100 * (v[2] - v[1]^2)^2
rosen_grad <- function(v) c(-2 * (1 - v[1]) - 400 * v[1] * (v[2] - v[1]^2),
                            200 * (v[2] - v[1]^2))
rosen_hess <- function(v) matrix(c(2 - 400 * v[2] + 1200 * v[1]^2, -400 * v[1],
                                   -400 * v[1], 200), 2)

test_that("external functions evaluate and differentiate through the chain rule", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 2)
  f <- opt_function(rosen, 2 * x, gradient = rosen_grad, hessian = rosen_hess,
                    name = "rosen")
  opt_objective(m, f + x[1])
  at <- c(0.3, -0.2)
  u <- 2 * at
  expect_equal(unname(opt_gradient(m, at)), 2 * rosen_grad(u) + c(1, 0))
  expect_equal(unname(as.matrix(opt_hessian(m, at))), 4 * rosen_hess(u))
  expect_match(format(f), "rosen")
})

test_that("supplied and finite difference gradients solve with labels", {
  skip_if_not_installed("nloptr")
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 2, start = c(-1.2, 1))
  opt_objective(m, opt_function(rosen, x, gradient = rosen_grad))
  r <- opt_solve(m)
  expect_equal(opt_value(r, x), c(1, 1), tolerance = 1e-5)
  expect_length(r$metadata$finite_difference, 0)
  opt_objective(m, opt_function(rosen, x, gradient = "finite_difference",
                                name = "rosen_fd"))
  r2 <- opt_solve(m)
  expect_equal(opt_value(r2, x), c(1, 1), tolerance = 1e-3)
  expect_identical(r2$metadata$finite_difference, "rosen_fd")
  expect_true(any(grepl("finite-difference", r2$warnings)))
})

test_that("bad external functions fail clearly", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 2)
  f <- opt_function(function(v) v, x, gradient = function(v) 1)
  opt_objective(m, f)
  expect_error(opt_gradient(m, c(1, 2)), "Gradient of 'f' has length 1")
  expect_error(opt_hessian(m, c(1, 2)), "No Hessian")
  expect_error(opt_function(rosen, x, gradient = "nope"))
})

test_that("jaxR gradients match the analytic gradient", {
  skip_if_not_installed("jaxR")
  skip_if(inherits(try(jaxR::devices(), silent = TRUE), "try-error"),
          "no XLA runtime for jaxR")
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 2)
  # jaxR indexing is 0-based, as in JAX: v[0] is the first element
  jf <- function(v) (1 - v[0])^2 + 100 * (v[1] - v[0]^2)^2
  opt_objective(m, opt_function(jf, x, gradient = "jaxR"))
  at <- c(0.3, -0.2)
  expect_equal(unname(opt_gradient(m, at)), rosen_grad(at), tolerance = 1e-6)
  expect_equal(optcompR:::eval_expr(m$objective, at), rosen(at), tolerance = 1e-12)
})
