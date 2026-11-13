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

test_that("Clarabel solves quadratic objectives with cross terms", {
  skip_if_not_installed("clarabel")
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 0)
  y <- opt_variable(m, "y", lower = 0)
  opt_constraint(m, x + y <= 1.5)
  opt_objective(m, x^2 + x * y + 2 * y^2 - 3 * x - 4 * y + 1)
  r <- opt_solve(m, solver = "clarabel")
  expect_equal(unname(opt_value(r, c(x, y))), c(0.875, 0.625), tolerance = 1e-6)
})

test_that("Clarabel keeps cone multipliers without ordinary rows", {
  skip_if_not_installed("clarabel")
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 2)
  t <- opt_variable(m, "t")
  opt_soc_constraint(m, x - 1, t)
  opt_objective(m, t)
  r <- opt_solve(m, solver = "clarabel")
  expect_length(r$extra$cone_duals, 3L)
})

test_that("ROI handles equality rows", {
  skip_if_not_installed("ROI")
  skip_if_not_installed("ROI.plugin.glpk")
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 2, lower = 0)
  opt_constraint(m, x[1] + x[2] == 2)
  opt_objective(m, x[1] + 2 * x[2])
  r <- opt_solve(m, solver = "roi", options = list(roi_solver = "glpk"))
  expect_equal(opt_objective_value(r), 2)
})

test_that("matrix indices select coordinates", {
  m <- opt_model()
  X <- opt_variable(m, "X", dim = c(2, 3))
  e <- X[cbind(c(1, 2), c(3, 1))]
  expect_equal(e$len, 2L)
  opt_objective(m, sum(e))
  g <- as.vector(opt_gradient(m, rep(0, 6)))
  expect_equal(which(g != 0), c(2L, 5L))
})

test_that("tiny indefinite quadratics are not convex", {
  m <- opt_model()
  x <- opt_variable(m, "x", lower = -1e6, upper = 1e6)
  y <- opt_variable(m, "y", lower = -1e6, upper = 1e6)
  opt_objective(m, 1e-14 * (x^2 + 4 * x * y + y^2))
  expect_identical(opt_classify(m)$convexity, "nonconvex")
})

test_that("quadratic rows without finite bounds do not affect convexity", {
  m <- opt_model()
  x <- opt_variable(m, "x")
  opt_constraint(m, x^2, lower = -Inf, upper = Inf)
  opt_objective(m, x)
  expect_identical(opt_classify(m)$convexity, "convex")
})

test_that("local solvers start inside the bounds", {
  skip_if_not_installed("nloptr")
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 0, upper = 10)
  opt_objective(m, x - log(x))
  r <- opt_solve(m, solver = "nloptr", time_limit = 20)
  expect_equal(opt_value(r, x), 1, tolerance = 1e-5)
})

test_that("a start where the objective is undefined is reported", {
  skip_if_not_installed("nloptr")
  m <- opt_model()
  x <- opt_variable(m, "x", lower = -1, upper = 1, start = 0)
  opt_objective(m, 1 / x^2 + x)
  expect_error(opt_solve(m, solver = "nloptr"), "opt_set_start")
})

test_that("nlpbb does not claim infeasibility after lost nodes", {
  skip_if_not_installed("nloptr")
  m <- opt_model()
  x <- opt_variable(m, "x", type = "integer", lower = -5, upper = 5)
  y <- opt_variable(m, "y")
  opt_objective(m, (x - 0.5)^2 - y)
  r <- suppressWarnings(opt_solve(m, solver = "nlpbb", time_limit = 30))
  expect_false(opt_status(r) %in% c("infeasible", "locally_infeasible"))
})

test_that("nlpbb incumbents satisfy the constraints after rounding", {
  skip_if_not_installed("nloptr")
  m <- opt_model()
  x <- opt_variable(m, "x", type = "integer", lower = 0, upper = 3)
  y <- opt_variable(m, "y", lower = -1e7, upper = 1e7)
  opt_constraint(m, y == 2e6 * x)
  opt_objective(m, (y - 1999999)^2 / 1e12 + exp(0 * y))
  r <- suppressWarnings(opt_solve(m, solver = "nlpbb", time_limit = 30))
  if (opt_has_solution(r)) {
    expect_true(opt_check(r, tol = 1e-6)$feasible)
  } else {
    succeed()
  }
})

test_that("solution derivatives refuse kinks and infeasible points", {
  m <- opt_model()
  x <- opt_variable(m, "x", start = 0)
  p <- opt_parameter(m, "p", value = 0)
  opt_objective(m, (x - p)^2 + abs(x))
  r <- opt_solve(m, solver = "nloptr")
  expect_error(opt_solution_jacobian(r, p), "abs")
})

test_that("near-bound interior points keep their derivatives", {
  skip_if_not_installed("clarabel")
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 0)
  p <- opt_parameter(m, "p", value = 1e-4)
  opt_objective(m, 1e4 * (x - p)^2)
  r <- opt_solve(m, solver = "clarabel")
  d <- opt_solution_jacobian(r, p)
  expect_equal(as.vector(d$solution), 1, tolerance = 1e-6)
})

test_that("indicator switches must stay binary", {
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 1, upper = 10)
  z <- opt_variable(m, "z", type = "binary")
  opt_indicator(m, z, x <= 0, active = 0)
  opt_objective(m, z)
  opt_set_type(m, z, "continuous")
  expect_error(opt_solve(m, solver = "highs"), "binary switch")
})

test_that("fixed relaxation refuses a changed model", {
  m <- opt_model()
  x <- opt_variable(m, "x", type = "integer", lower = 0, upper = 5)
  p <- opt_parameter(m, "p", value = 1)
  opt_constraint(m, x >= p)
  opt_objective(m, x)
  r <- opt_solve(m, solver = "highs")
  opt_set_parameter(m, p, 3)
  expect_error(opt_fixed_relaxation(r), "changed")
})

test_that("automatic selection skips solvers that cannot reformulate", {
  skip_if_not_installed("scip")
  m <- opt_model("max")
  x <- opt_variable(m, "x", lower = 0)
  y <- opt_variable(m, "y", lower = 0)
  opt_constraint(m, x + y <= 4)
  opt_complementarity(m, x, y)
  opt_objective(m, x + 2 * y)
  r <- opt_solve(m)
  expect_identical(r$solver, "scip")
})

test_that("ignored controls are reported", {
  skip_if_not_installed("osqp")
  m <- opt_model()
  x <- opt_variable(m, "x", start = 1)
  opt_objective(m, (x - 2)^2)
  r <- opt_solve(m, solver = "osqp")
  expect_true(any(grepl("start values", r$warnings)))
})

test_that("element names are unique across blocks", {
  m <- opt_model()
  opt_variable(m, "x", dim = 2)
  expect_error(opt_variable(m, "x[1]"), "already used")
})

test_that("mean of an expression", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 4)
  opt_objective(m, mean(x))
  expect_equal(as.vector(opt_gradient(m, rep(0, 4))), rep(0.25, 4))
})
