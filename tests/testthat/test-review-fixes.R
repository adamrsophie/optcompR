test_that("binary variables keep 0/1 bounds after opt_set_bounds", {
  m <- opt_model("max")
  x <- opt_variable(m, "x", type = "binary")
  opt_set_bounds(m, x, upper = 5)
  opt_objective(m, x)
  expect_equal(opt_value(opt_solve(m, solver = "highs"), x), 1)
})

test_that("slightly negative curvature is not called convex", {
  m <- opt_model()
  x <- opt_variable(m, "x", lower = -1e6, upper = 1e6)
  y <- opt_variable(m, "y", lower = -1e6, upper = 1e6)
  opt_objective(m, (x - y)^2 - 1e-10 * (x + y)^2)
  expect_false(identical(opt_classify(m)$convexity, "convex"))
  expect_error(opt_solve(m, solver = "highs"), "not proven convex|nonconvex")
})

test_that("NLopt stopval is an objective limit, not optimal", {
  expect_equal(optcompR:::nloptr_termination(2L), "objective_limit")
  expect_equal(optcompR:::nloptr_termination(3L), "optimal")
})

test_that("impossible infinite row bounds are rejected", {
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 0, upper = 1)
  expect_error(opt_constraint(m, x, lower = Inf, upper = Inf), "never hold")
  expect_error(opt_constraint(m, x, lower = -Inf, upper = -Inf), "never hold")
})

test_that("log() honors base", {
  m <- opt_model()
  x <- opt_variable(m, "x")
  opt_objective(m, log(x, base = 10))
  expect_equal(unname(opt_gradient(m, 100)), 1 / (100 * log(10)))
  expect_error(log(x, base = -1), "base")
})

test_that("nlpbb proves a convex MIQP without NaN pruning", {
  skip_if_not_installed("nloptr")
  m <- opt_model()
  x <- opt_variable(m, "x", type = "integer", lower = 0, upper = 2)
  opt_objective(m, (x - 0.4)^2)
  r <- opt_solve(m, solver = "nlpbb")
  expect_equal(opt_value(r, x), 0, tolerance = 1e-6)
  expect_equal(opt_status(r), "optimal")
})

test_that("SOS2 with two nonadjacent nonzeros is a violation", {
  sf <- list(cones = list(), indicators = list(),
             sos = list(list(name = "s", type = 2L, cols = 1:3)))
  v <- optcompR:::special_violations(sf, c(1, 0, 1), 1e-9)
  expect_gt(v$violation, 0.5)
  v2 <- optcompR:::special_violations(sf, c(0, 1, 1), 1e-9)
  expect_equal(v2$violation, 0)
})

test_that("vector powers with exponents 0 and 1 have finite derivatives", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 2)
  opt_objective(m, sum(x^c(0, 1)) + sum(x^c(2, 2)))
  expect_equal(as.vector(opt_gradient(m, c(0, 0))), c(0, 1))
  expect_equal(as.matrix(opt_hessian(m, c(0, 0))), diag(2, 2),
               ignore_attr = TRUE)
})

test_that("row times column products of matrix expressions work", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 2)
  e <- t(x) %*% (diag(2) %*% x)
  opt_objective(m, e)
  expect_equal(as.vector(opt_gradient(m, c(1, 2))), c(2, 4))
})

test_that("solution derivatives need a KKT point", {
  skip_if_not_installed("nloptr")
  m <- opt_model()
  x <- opt_variable(m, "x", start = 0)
  p <- opt_parameter(m, "p", value = 1)
  opt_objective(m, (x - p)^2)
  r <- suppressWarnings(opt_solve(m, solver = "nloptr",
                                  options = list(maxeval = 1L)))
  expect_error(opt_solution_jacobian(r, p), "optimum|KKT")
})

test_that("bound-active solutions have zero solution derivatives", {
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 0)
  p <- opt_parameter(m, "p", value = 1)
  opt_objective(m, p * x)
  r <- opt_solve(m, solver = "highs")
  d <- opt_solution_jacobian(r, p)
  expect_equal(as.vector(d$solution), 0)
  expect_equal(unname(d$objective), 0)
})

test_that("indefinite KKT points are refused", {
  m <- opt_model()
  x <- opt_variable(m, "x", start = 0)
  p <- opt_parameter(m, "p", value = 0)
  opt_objective(m, -(x - p)^2)
  r <- opt_solve(m, solver = "nloptr")
  if (opt_has_solution(r) && abs(opt_value(r, x)) < 1e-8) {
    expect_error(opt_solution_jacobian(r, p), "second-order")
  } else {
    succeed()
  }
})

test_that("estimated multipliers respect sign restrictions", {
  m <- opt_model()
  x <- opt_variable(m, "x")
  opt_constraint(m, x <= 0, name = "a")
  opt_constraint(m, x >= 0, name = "b")
  opt_objective(m, x)
  r <- opt_solve(m, solver = "highs")
  r$row_dual <- NULL
  r$col_dual <- NULL
  expect_true(opt_check(r)$kkt$pass)
})

test_that("fixed relaxation resolves indicators", {
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 0, upper = 8)
  on <- opt_variable(m, "on", type = "binary")
  opt_indicator(m, on, x <= 0, active = 0)
  cap <- opt_constraint(m, x >= 3, name = "need")
  opt_objective(m, 2 * x + 10 * on)
  r <- suppressWarnings(opt_solve(m, solver = "highs"))
  fr <- opt_fixed_relaxation(r)
  expect_false(fr$classification$integer)
  expect_equal(unname(opt_dual(fr, cap)), 2)
  expect_length(m$special$indicator, 1L)
})

test_that("parameterized affine expressions enter cones", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 2)
  t <- opt_variable(m, "t")
  p <- opt_parameter(m, "p", value = 2)
  expect_silent(opt_soc_constraint(m, p * x, t))
})

test_that("constraint references are tied to their model", {
  a <- opt_model()
  xa <- opt_variable(a, "x")
  ca <- opt_constraint(a, xa <= 1, name = "c1")
  b <- opt_model()
  xb <- opt_variable(b, "x")
  opt_constraint(b, xb <= 1, name = "c1")
  expect_error(opt_delete_constraint(b, ca), "different model")
})

test_that("scenario values are validated", {
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 0)
  d <- opt_parameter(m, "d", value = 1)
  opt_constraint(m, x >= d)
  opt_objective(m, x)
  expect_error(opt_scenarios(m, d, list(NA_real_)))
})

test_that("subnormal coefficients convert to exact rationals", {
  txt <- optcompR:::z3_number_text(2^-1074, "binary")
  expect_match(txt, "^1/[0-9]+$")
  expect_false(grepl("Inf", txt))
  expect_equal(nchar(sub("^1/", "", txt)), 324L)
})

test_that("repeated robust and CVaR builders get distinct names", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 2, lower = 0, upper = 1)
  opt_robust_constraint(m, x, a0 = c(1, 1), b = 3, delta = 0.1)
  expect_silent(opt_robust_constraint(m, x, a0 = c(1, 2), b = 4, delta = 0.1))
  P <- diag(2)
  opt_robust_constraint(m, x, a0 = c(1, 1), b = 3, uncertainty = "ellipsoid",
                        P = P)
  expect_silent(opt_robust_constraint(m, x, a0 = c(1, 1), b = 5,
                                      uncertainty = "ellipsoid", P = P))
  L <- opt_cvar(m, x, alpha = 0.5)
  expect_silent(opt_cvar(m, x, alpha = 0.9))
})

test_that("a variable far from a distant opposite bound is not active", {
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 0, upper = 1e6)
  p <- opt_parameter(m, "p", value = 1)
  opt_objective(m, (x - p)^2)
  r <- opt_solve(m, solver = "highs")
  d <- opt_solution_jacobian(r, p)
  expect_equal(as.vector(d$solution), 1, tolerance = 1e-6)
})

test_that("fixed relaxation keeps cone constraints", {
  skip_if_not_installed("scip")
  m <- opt_model("max")
  x <- opt_variable(m, "x", dim = 2, lower = -5, upper = 5)
  t <- opt_variable(m, "t", lower = 0, upper = 1)
  z <- opt_variable(m, "z", type = "binary")
  opt_soc_constraint(m, x, t)
  opt_indicator(m, z, t <= 0.5, active = 0)
  opt_objective(m, sum(x) - 0.1 * z)
  r <- suppressWarnings(opt_solve(m, solver = "scip"))
  fr <- suppressWarnings(opt_fixed_relaxation(r, solver = "scip"))
  expect_lte(sqrt(sum(opt_value(fr, x)^2)), opt_value(fr, t) + 1e-6)
})
