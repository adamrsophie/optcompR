test_that("symbolic expressions convert exactly", {
  skip_if_not_installed("symcompR")
  S <- symcompR::Symbol
  m <- opt_model()
  x <- opt_variable(m, "x")
  y <- opt_variable(m, "y")
  s <- (S("x") - 1)^2 + symcompR::Rational(1, 3) * S("x") * S("y") +
    exp(S("y")) / 2 + sqrt(S("x"))
  e <- opt_from_sym(m, s)
  f <- function(x, y) (x - 1)^2 + x * y / 3 + exp(y) / 2 + sqrt(x)
  z <- c(1.7, -0.4)
  expect_equal(optcompR:::eval_expr(e, z), f(1.7, -0.4))
  # symbolic derivative agrees with optcompR's exact gradient
  opt_objective(m, e)
  dsx <- opt_from_sym(m, diff(s, S("x")))
  expect_equal(optcompR:::eval_expr(dsx, z), unname(opt_gradient(m, z)[1]))
})

test_that("quadratic symbolic objectives classify as QP and solve", {
  skip_if_not_installed("symcompR")
  S <- symcompR::Symbol
  m <- opt_model()
  x <- opt_variable(m, "x")
  y <- opt_variable(m, "y")
  opt_objective(m, opt_from_sym(m, (S("x") - 1)^2 + (S("y") - 2)^2))
  opt_constraint(m, opt_from_sym(m, S("x") + S("y")) >= 4)
  expect_identical(opt_classify(m)$class, "QP")
  expect_equal(opt_objective_value(opt_solve(m)), 0.5, tolerance = 1e-6)
})

test_that("round trips through symcompR preserve the function", {
  skip_if_not_installed("symcompR")
  m <- opt_model()
  x <- opt_variable(m, "x")
  y <- opt_variable(m, "y")
  p <- opt_parameter(m, "p", value = 2)
  e <- 3 * x^2 - x * y + p * exp(y) / (1 + x^2) + 0.5
  back <- opt_from_sym(m, opt_to_sym(e))
  z <- c(0.3, 1.1, 2)
  expect_equal(optcompR:::eval_expr(back, z), optcompR:::eval_expr(e, z))
  expect_error(opt_from_sym(m, symcompR::Symbol("w")), "map")
  w <- opt_from_sym(m, symcompR::Symbol("w") + 1, map = list(w = x))
  expect_equal(optcompR:::eval_expr(w, z), 1.3)
})
