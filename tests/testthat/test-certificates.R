# Bounds, optimality claims and derivatives checked against known answers
# on badly scaled and edge-case models.

test_that("branch and bound keeps a node whose cut LP drops coefficients", {
  skip_if_not_installed("nloptr")
  skip_if_not_installed("highs")
  m <- opt_model()
  x <- opt_variable(m, "x", type = "integer", lower = 0, upper = 1e7,
                    start = 5e6)
  opt_constraint(m, 1e-5 * exp(-x / 1e7) <= 1e-5 * exp(-0.5))
  opt_objective(m, (x / 1e7 - 0.5)^2)
  r <- suppressWarnings(opt_solve(m, solver = "nlpbb"))
  expect_true(opt_has_solution(r))
  expect_lt(opt_objective_value(r), 1e-9)
})

test_that("open nodes report certified bounds in both senses", {
  skip_if_not_installed("nloptr")
  skip_if_not_installed("highs")
  for (sense in c("min", "max")) {
    m <- opt_model(sense)
    x <- opt_variable(m, "x", type = "integer", lower = 0, upper = 1e7,
                      start = 0.5)
    opt_constraint(m, x^2 <= 1e14)
    sg <- if (sense == "min") 1 else -1
    opt_objective(m, sg * (x / 1e7 - 0.6)^2)
    r <- opt_solve(m, solver = "nlpbb", options = list(max_nodes = 1L))
    b <- opt_objective_bound(r)
    # x = 6e6 attains 0
    if (sense == "min") expect_lte(b, 1e-9) else expect_gte(b, -1e-9)
  }
})

test_that("a nearly integral relaxation does not close its node", {
  skip_if_not_installed("nloptr")
  skip_if_not_installed("highs")
  m <- opt_model()
  x <- opt_variable(m, "x", type = "binary")
  y <- opt_variable(m, "y")
  opt_constraint(m, y >= 1 - 2e6 * x)
  opt_constraint(m, y >= 0.1 * x)
  opt_objective(m, 20000 + y + 0.01 * x^2)
  r <- opt_solve(m, solver = "nlpbb")
  expect_equal(opt_objective_value(r), 20000.11, tolerance = 1e-8)
  expect_lte(opt_objective_bound(r), 20000.11 + 1e-6)
})

test_that("one-sided boxes keep finite bounds; heuristic nodes close", {
  skip_if_not_installed("nloptr")
  skip_if_not_installed("highs")
  m <- opt_model()
  x <- opt_variable(m, "x", type = "integer", lower = 0, upper = 2)
  y <- opt_variable(m, "y", lower = 0)
  opt_objective(m, (x - 0.4)^2 + y)
  expect_true(is.finite(oa_lower_bound(opt_standard_form(m), c(0.4, 0), 1)$bound))
  r <- opt_solve(m, solver = "nlpbb")
  expect_equal(opt_status(r), "optimal")
  expect_equal(opt_objective_value(r), 0.16, tolerance = 1e-8)

  m <- opt_model()
  x <- opt_variable(m, "x", type = "integer", lower = 0)
  opt_objective(m, (x + 1)^2)
  r <- opt_solve(m, solver = "nlpbb", options = list(max_nodes = 3L))
  expect_equal(unname(r$x), 0)

  m <- opt_model()
  x <- opt_variable(m, "x", type = "integer", lower = -100, upper = 100,
                    start = 1)
  opt_objective(m, (x^2 - 1)^2)
  r <- opt_solve(m, solver = "nlpbb", options = list(max_nodes = 3L))
  expect_equal(opt_status(r), "locally_optimal")
  expect_identical(opt_solver_capabilities("nlpbb")$global, character())
})

test_that("rank deficient PSD rows stay convex", {
  set.seed(1)
  for (i in 1:50) {
    a <- runif(3, -3, 3)
    m <- opt_model()
    x <- opt_variable(m, "x", dim = 3)
    opt_constraint(m, (a[1] * x[1] + a[2] * x[2] + a[3] * x[3])^2 <= 1)
    opt_constraint(m, (a[1] * x[1] + a[2] * x[2])^2, lower = -1, upper = 1)
    opt_objective(m, (a[1] * x[1] + a[2] * x[2] - 1)^2)
    expect_equal(opt_classify(m)$convexity, "convex")
  }
})

test_that("tiny negative curvature and fixed variables are not ignored", {
  m <- opt_model()
  x <- opt_variable(m, "x", lower = -2, upper = 2)
  y <- opt_variable(m, "y", lower = 1e7, upper = 1e7)
  opt_constraint(m, x^2 - 1e-16 * y^2, lower = 0, upper = 4)
  opt_objective(m, x)
  expect_false(opt_classify(m)$convexity == "convex")

  m <- opt_model()
  x <- opt_variable(m, "x", lower = -2, upper = 2)
  y <- opt_variable(m, "y", lower = -1e14, upper = -1e14)
  opt_constraint(m, x^2 + 1e-14 * y, lower = 0, upper = 1)
  opt_objective(m, x)
  expect_equal(opt_classify(m)$convexity, "nonconvex")

  # the box removes the negative branch: no witness of nonconvexity
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 0, upper = 2)
  opt_constraint(m, x^2, lower = 1, upper = 4)
  opt_objective(m, x)
  expect_false(opt_classify(m)$convexity == "nonconvex")
})

test_that("an NLopt tolerance stop is not a global optimum", {
  skip_if_not_installed("nloptr")
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 0, upper = 1e7, start = 0.5)
  opt_objective(m, (x / 1e7 - 0.6)^2)
  r <- opt_solve(m, solver = "nloptr")
  expect_false(opt_status(r) == "optimal")
})

test_that("backends must return a whole point with a feasible claim", {
  opt_register_solver("cert_missing", function(sf, cls, control)
    list(termination = "optimal", primal_status = "feasible_point",
         x = NA_real_), classes = "LP", global = "LP", version = "test")
  on.exit(opt_unregister_solver("cert_missing"), add = TRUE)
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 0, upper = 1)
  opt_objective(m, x)
  expect_error(opt_solve(m, solver = "cert_missing"), "complete finite")
})

test_that("opt_check reports a point where a row is undefined", {
  opt_register_solver("cert_bad_row", function(sf, cls, control)
    list(termination = "optimal", primal_status = "feasible_point", x = 0),
    classes = "NLP", version = "test")
  on.exit(opt_unregister_solver("cert_bad_row"), add = TRUE)
  m <- opt_model()
  x <- opt_variable(m, "x")
  opt_constraint(m, sqrt(x - 1) <= 2)
  opt_objective(m, x)
  r <- suppressWarnings(opt_solve(m, solver = "cert_bad_row"))
  d <- suppressWarnings(opt_check(r))
  expect_false(d$feasible)
  expect_false(d$kkt$pass)
})

test_that("validation uses the strictest requested accuracy", {
  e <- list(feasibility_tol = 1e-2, accuracy_options = c("eps_abs", "eps_rel"))
  v <- validation_policy(e, list(options = list(eps_abs = 1e-9, eps_rel = 1e-3)))
  expect_equal(v$tol, 1e-7)
})

test_that("infinite bounds recycle and print", {
  m <- opt_model()
  x <- opt_variable(m, "x")
  expect_output(print(x + 1 <= Inf), "x <= Inf")
  opt_constraint(m, x <= c(Inf, 2))
  opt_objective(m, -x)
  skip_if_not_installed("highs")
  expect_equal(unname(opt_solve(m, solver = "highs")$x), 2)
  expect_error(x <= c(NA, Inf), "missing")
})

test_that("constant implications hold or are rejected", {
  m <- opt_model()
  x <- opt_variable(m, "x", type = "binary")
  expect_null(opt_implies(m, 0, 1))
  expect_error(opt_implies(m, 1, 0), "never hold")
})

test_that("ranged rows report their active side; bound-only LPs range", {
  skip_if_not_installed("highs")
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 0, upper = 10)
  opt_constraint(m, x, lower = 2, upper = 6, name = "range")
  opt_objective(m, x + 10)
  rw <- opt_sensitivity(opt_solve(m, solver = "highs"))$rows
  expect_equal(rw$side, "lower")
  expect_equal(c(rw$bound_lower, rw$bound_upper), c(0, 6))

  m <- opt_model()
  x <- opt_variable(m, "x", lower = 0, upper = 1)
  opt_objective(m, x)
  s <- opt_sensitivity(opt_solve(m, solver = "highs"))
  expect_equal(c(s$costs$lower, s$costs$upper), c(0, Inf))

  m <- opt_model("max")
  x <- opt_variable(m, "x", lower = 1, upper = 2)
  opt_constraint(m, 0 * x <= 1)
  opt_objective(m, 3 * x)
  s <- opt_sensitivity(opt_solve(m, solver = "highs"))
  expect_equal(c(s$costs$lower, s$costs$upper), c(0, Inf))
  expect_equal(s$rows$bound_upper, Inf)
})

test_that("a rebuilt fixed relaxation ignores later edits", {
  skip_if_not_installed("highs")
  m <- opt_model()
  x <- opt_variable(m, "x", type = "integer", lower = 0, upper = 3)
  y <- opt_variable(m, "y", lower = 0)
  p <- opt_parameter(m, "p", 1)
  opt_constraint(m, y >= p)
  opt_objective(m, sqrt(p) * x + 2 * y)
  r <- opt_solve(m, solver = "highs")
  w <- opt_variable(m, "w", type = "integer", lower = 4, upper = 6)
  f <- opt_fixed_relaxation(r, solver = "highs")
  expect_error(opt_reduced_cost(f, w), "added after")
  expect_equal(f$classification$class, "LP")
  expect_equal(names(f$x), c("x", "y"))
  opt_set_parameter(m, p, -1)
  f <- opt_fixed_relaxation(r, solver = "highs")
  expect_equal(opt_value(f, y), 1)
  expect_equal(unname(opt_reduced_cost(f, y)), 0)
  expect_true(is.list(opt_solution_jacobian(f, p)))
})
