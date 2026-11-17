test_that("parameter exponents stay power nodes until values are known", {
  m <- opt_model()
  x <- opt_variable(m, "x", lower = -5, upper = 5)
  p <- opt_parameter(m, "p", value = 2)
  cc <- opt_parameter(m, "c", value = 3)
  e <- x^p
  expect_equal(e$op, "powp")
  opt_objective(m, e + cc * x)
  expect_equal(opt_classify(m)$class, "QP")
  r <- opt_solve(m, solver = "nloptr")
  expect_equal(unname(opt_value(r, x)), -1.5, tolerance = 1e-5)
  expect_equal(unname(opt_solution_jacobian(r, cc)$solution[1, 1]), -0.5,
               tolerance = 1e-5)
  expect_error(opt_solution_jacobian(r, p), "not differentiable in p")
  opt_set_parameter(m, p, 4)
  expect_equal(opt_classify(m)$class, "NLP")
})

test_that("derivatives of a^p in the exponent are exact for a positive base", {
  m <- opt_model()
  y <- opt_variable(m, "y", lower = 0.5, upper = 3)
  q <- opt_parameter(m, "q", value = 1.5)
  opt_objective(m, (y^q - 2)^2)
  r <- opt_solve(m, solver = "nloptr")
  y0 <- 2^(1 / 1.5)
  expect_equal(unname(opt_value(r, y)), y0, tolerance = 1e-5)
  expect_equal(unname(opt_solution_jacobian(r, q)$solution[1, 1]),
               -log(y0) / 1.5 * y0, tolerance = 1e-4)
})

test_that("parameter powers round trip through symcompR", {
  skip_if_not_installed("symcompR")
  m <- opt_model()
  x <- opt_variable(m, "x")
  p <- opt_parameter(m, "p", value = 3)
  s <- opt_to_sym(x^p)
  back <- opt_from_sym(m, s)
  expect_equal(back$op, "powp")
  expect_equal(format(back), format(x^p))
})
