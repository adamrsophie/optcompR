# Global optimality claims must survive a better feasible witness, and
# reported bounds must never exceed a known feasible objective.

test_that("bilinear two-sided rows are never proved convex", {
  for (perm in list(c(1, 2), c(2, 1))) {
    m <- opt_model()
    x <- opt_variable(m, "x", dim = 3, lower = -2, upper = 2,
                      start = c(.5, .5, .5))
    opt_constraint(m, x[perm[1]] * x[perm[2]], lower = .1, upper = 1)
    opt_objective(m, sum(x[1:2]))
    expect_false(opt_classify(m)$convexity == "convex")
  }
  skip_if_not_installed("nloptr")
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 2, lower = -2, upper = 2, start = c(.5, .5))
  opt_constraint(m, x[1] * x[2], lower = .1, upper = 1)
  opt_objective(m, sum(x))
  r <- opt_solve(m, solver = "nloptr")
  # (-2, -0.5) is feasible with objective -2.5
  expect_false(opt_status(r) == "optimal" && opt_objective_value(r) > -2.5 + 1e-5)
})

test_that("a positive lower bound on a square is never treated as redundant", {
  for (sc in c(1, 1e-3, 1e4)) {
    m <- opt_model()
    x <- opt_variable(m, "x", lower = -2, upper = 2, start = .5)
    opt_constraint(m, sc * x^2, lower = sc * 1e-5, upper = sc)
    opt_objective(m, x)
    expect_false(opt_classify(m)$convexity == "convex")
    # shifted square
    m <- opt_model()
    x <- opt_variable(m, "x", lower = -2, upper = 4)
    opt_constraint(m, sc * (x - 1)^2, lower = sc * 1e-6, upper = sc)
    opt_objective(m, x)
    expect_false(opt_classify(m)$convexity == "convex")
  }
  # exact redundancy is still recognized
  m <- opt_model()
  x <- opt_variable(m, "x", lower = -2, upper = 2)
  opt_constraint(m, x^2, lower = 0, upper = 1)
  opt_objective(m, x)
  expect_identical(opt_classify(m)$convexity, "convex")
  skip_if_not_installed("nloptr")
  m <- opt_model()
  x <- opt_variable(m, "x", lower = -2, upper = 2, start = .5)
  opt_constraint(m, x^2, lower = 1e-5, upper = 1)
  opt_objective(m, x)
  r <- opt_solve(m, solver = "nloptr")
  expect_false(opt_status(r) == "optimal" && opt_objective_value(r) > -1 + 1e-5)
})

test_that("NLP branch and bound bounds are certified, not primal values", {
  skip_if_not_installed("nloptr")
  m <- opt_model()
  x <- opt_variable(m, "x", type = "integer", lower = -1e7, upper = 1e7)
  opt_objective(m, (x / 1e7 - .6)^2)
  r <- opt_solve(m, solver = "nlpbb", mip_gap = 0)
  b <- opt_objective_bound(r)
  # x = 6e6 gives 0
  expect_true(is.na(b) || b <= 1e-8)
  expect_false(opt_status(r) == "optimal" && opt_objective_value(r) > 1e-8)
  expect_false(isTRUE(opt_gap(r) == 0) && opt_objective_value(r) > 1e-8)
  # maximization: the bound is an upper bound
  m <- opt_model("max")
  x <- opt_variable(m, "x", type = "integer", lower = -1e7, upper = 1e7)
  opt_objective(m, -(x / 1e7 - .6)^2)
  r <- opt_solve(m, solver = "nlpbb", mip_gap = 0)
  b <- opt_objective_bound(r)
  expect_true(is.na(b) || b >= -1e-8)
  # an honest gap is kept when the requested gap stops the search early
  m <- opt_model()
  x <- opt_variable(m, "x", type = "integer", lower = 0, upper = 1)
  opt_objective(m, (x - .6)^2)
  r <- opt_solve(m, solver = "nlpbb", mip_gap = .5)
  expect_true(opt_objective_bound(r) <= 0.16 + 1e-9)
  r <- opt_solve(m, solver = "nlpbb")
  expect_identical(opt_status(r), "optimal")
  expect_equal(opt_objective_value(r), 0.16, tolerance = 1e-8)
})
