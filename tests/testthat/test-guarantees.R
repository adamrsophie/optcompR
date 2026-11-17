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
  # a parameter at zero stores a zero coefficient in the row
  m <- opt_model()
  x <- opt_variable(m, "x", lower = -2, upper = 2)
  y <- opt_variable(m, "y", lower = -2, upper = 2)
  p <- opt_parameter(m, "p", 0)
  opt_constraint(m, x^2 + p * x + y, lower = -1, upper = 1)
  opt_objective(m, x)
  expect_false(opt_classify(m)$convexity == "convex")
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

test_that("NLP branch and bound certifies through convex nonlinear rows", {
  skip_if_not_installed("nloptr")
  m <- opt_model()
  x <- opt_variable(m, "x", type = "integer", lower = 0, upper = 5)
  y <- opt_variable(m, "y", type = "integer", lower = 0, upper = 5)
  opt_constraint(m, x^2 + y^2 <= 8)
  opt_objective(m, (x - 2.6)^2 + (y - 2.6)^2)
  r <- opt_solve(m, solver = "nlpbb")
  expect_identical(opt_status(r), "optimal")
  expect_equal(unname(r$x), c(2, 2))
  # enumeration of the integer points confirms the optimum 0.72
  g <- expand.grid(x = 0:5, y = 0:5)
  g <- g[g$x^2 + g$y^2 <= 8, ]
  expect_equal(min((g$x - 2.6)^2 + (g$y - 2.6)^2), opt_objective_value(r),
               tolerance = 1e-8)
  expect_lte(opt_objective_bound(r), opt_objective_value(r) + 1e-12)
})

test_that("a local mixed-integer backend is not promoted by a convex relaxation", {
  name <- "test_local_miqp"
  opt_register_solver(name, function(sf, cls, control) {
    list(termination = "optimal", primal_status = "feasible_point",
         dual_status = "not_available", x = rep(0, length(sf$c)),
         objective = sf$offset)
  }, classes = "MIQP", global = character(), version = "test")
  on.exit(opt_unregister_solver(name), add = TRUE)
  m <- opt_model()
  x <- opt_variable(m, "x", type = "integer", lower = 0, upper = 1)
  opt_objective(m, (x - .6)^2)
  expect_identical(opt_status(opt_solve(m, solver = name)), "locally_optimal")
})

test_that("an undefined objective at the returned point is not a solution", {
  name <- "test_zero_nlp"
  opt_register_solver(name, function(sf, cls, control) {
    list(termination = "optimal", primal_status = "feasible_point",
         dual_status = "not_available", x = rep(0, length(sf$c)), objective = 0)
  }, classes = "NLP", version = "test")
  on.exit(opt_unregister_solver(name), add = TRUE)
  for (f in list(log, function(x) sqrt(x - 1))) {
    m <- opt_model()
    x <- opt_variable(m, "x", lower = 0)
    opt_objective(m, f(x))
    r <- opt_solve(m, solver = name)
    expect_false(opt_has_solution(r))
    expect_identical(opt_status(r), "numerical_error")
    expect_false(anyNA(r$x))
  }
})

test_that("the KKT verdict requires primal feasibility", {
  name <- "test_zero_qp"
  opt_register_solver(name, function(sf, cls, control) {
    list(termination = "optimal", primal_status = "feasible_point",
         dual_status = "feasible_point", x = rep(0, length(sf$c)), objective = 0,
         row_dual = rep(0, nrow(sf$A)), col_dual = rep(0, length(sf$c)))
  }, classes = "QP", version = "test")
  on.exit(opt_unregister_solver(name), add = TRUE)
  m <- opt_model()
  x <- opt_variable(m, "x")
  opt_constraint(m, x >= 1)
  opt_objective(m, x^2)
  d <- opt_check(opt_solve(m, solver = name))
  expect_false(d$feasible)
  expect_false(d$kkt$pass)
  expect_true(d$kkt$dual_pass)
  expect_false(d$kkt$primal_feasibility)
})

test_that("a saved OSQP model rebuilds its workspace", {
  skip_if_not_installed("osqp")
  m <- opt_model()
  x <- opt_variable(m, "x", lower = -5, upper = 5)
  p <- opt_parameter(m, "p", 1)
  opt_objective(m, (x - p)^2)
  opt_solve(m, solver = "osqp")
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f), add = TRUE)
  saveRDS(m, f)
  m2 <- readRDS(f)
  r <- opt_solve(m2, solver = "osqp")
  expect_equal(as.numeric(r$x), 1, tolerance = 1e-5)
  expect_false(isTRUE(r$extra$reused))
  r <- opt_solve(m2, solver = "osqp")
  expect_true(isTRUE(r$extra$reused))
  opt_set_parameter(m2, "p", 2)
  expect_equal(as.numeric(opt_solve(m2, solver = "osqp")$x), 2,
               tolerance = 1e-5)
})

test_that("a failed complementarity call leaves the model unchanged", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 2, lower = 0, upper = 1)
  y <- opt_variable(m, "y", dim = 2, lower = 0, upper = 1)
  opt_sos1(m, c(x[2], y[2]), name = "pair[2]")
  before <- opt_special_constraints(m)
  expect_error(opt_complementarity(m, x, y, name = "pair"))
  expect_identical(opt_special_constraints(m), before)
  opt_constraint(m, x[1] == 1)
  opt_constraint(m, y[1] == 1)
  opt_objective(m, 0 * x[1])
  expect_identical(opt_status(opt_solve(m, solver = "highs")), "optimal")
})

test_that("retained helper results form their own fixed relaxations", {
  m <- opt_model()
  x <- opt_variable(m, "x", type = "integer", lower = 0, upper = 3)
  opt_objective(m, x)
  r <- opt_weighted(m, list(x, 2 * x), c(1, 1), solver = "highs")
  fr <- opt_fixed_relaxation(r, solver = "highs")
  expect_equal(opt_objective_value(fr), opt_objective_value(r))
  # the user's handles index the rebuilt result
  m <- opt_model("max")
  x <- opt_variable(m, "x", type = "integer", lower = 0, upper = 3)
  y <- opt_variable(m, "y", lower = 0)
  cap <- opt_constraint(m, x + y <= 4.5)
  opt_objective(m, 3 * x + 2 * y)
  r <- opt_weighted(m, list(3 * x + 2 * y, x), c(1, 1), solver = "highs")
  fr <- opt_fixed_relaxation(r, solver = "highs")
  expect_equal(unname(opt_value(fr, x)), 3)
  expect_equal(unname(opt_value(fr, y)), 1.5)
  expect_equal(unname(opt_dual(fr, cap)), 2)
  m <- opt_model()
  x <- opt_variable(m, "x", type = "integer", lower = 0, upper = 10)
  p <- opt_parameter(m, "p", 2)
  opt_constraint(m, x >= p)
  opt_objective(m, x)
  before <- opt_standard_form(m)
  ss <- opt_scenarios(m, p, c(2, 3), keep = TRUE, solver = "highs")
  for (r in attr(ss, "results")) {
    fr <- opt_fixed_relaxation(r, solver = "highs")
    expect_equal(opt_objective_value(fr), opt_objective_value(r))
  }
  # the live model is untouched
  after <- opt_standard_form(m)
  drop <- c("key", "structure_key")
  expect_identical(before[setdiff(names(before), drop)],
                   after[setdiff(names(after), drop)])
})

test_that("LP row ranges agree with re-solves at their endpoints", {
  rhs0 <- c(capacity = 4, xmax = 2, ymax = 3, low = 1)
  build <- function(rhs) {
    m <- opt_model("max")
    x <- opt_variable(m, "x", lower = 0)
    y <- opt_variable(m, "y", lower = 0)
    opt_constraint(m, x + y <= rhs[["capacity"]], name = "capacity")
    opt_constraint(m, x <= rhs[["xmax"]], name = "xmax")
    opt_constraint(m, y <= rhs[["ymax"]], name = "ymax")
    opt_constraint(m, x + 3 * y >= rhs[["low"]], name = "low")
    opt_objective(m, 3 * x + 2 * y)
    m
  }
  r <- opt_solve(build(rhs0), solver = "highs")
  rows <- opt_sensitivity(r)$rows
  for (i in seq_len(nrow(rows))) {
    for (end in c("lower", "upper")) {
      b <- rows[[paste0("bound_", end)]][i]
      if (!is.finite(b)) next
      rhs <- rhs0
      rhs[[rows$row[i]]] <- b
      actual <- opt_objective_value(opt_solve(build(rhs), solver = "highs"))
      expect_equal(rows[[paste0("objective_at_", end)]][i], actual,
                   tolerance = 1e-7, info = paste(rows$row[i], end))
    }
  }
  z <- rows[rows$row == "ymax", ]
  expect_identical(z$side, "upper")
  expect_equal(z$bound_lower, 2)
  expect_identical(z$bound_upper, Inf)
})

test_that("a strict accuracy request tightens result validation", {
  skip_if_not_installed("scs")
  m <- opt_model()
  x <- opt_variable(m, "x")
  t <- opt_variable(m, "t")
  opt_constraint(m, x == 1)
  opt_soc_constraint(m, x, t)
  opt_objective(m, t)
  r <- opt_solve(m, solver = "scs",
                 options = list(max_iters = 15L, eps_abs = 1e-9, eps_rel = 1e-9))
  expect_identical(r$validation$source, "requested solver accuracy")
  expect_equal(r$validation$tol, 1e-7)
  if (!opt_check(r, tol = 1e-6)$feasible) {
    expect_false(opt_has_solution(r))
  }
  # a strict request that converges still gives a solution
  r <- opt_solve(m, solver = "scs", options = list(eps_abs = 1e-9, eps_rel = 1e-9))
  expect_identical(opt_status(r), "optimal")
  expect_true(opt_has_solution(r))
  r <- opt_solve(m, solver = "scs")
  expect_identical(r$validation$source, "solver default")
  expect_true(opt_has_solution(r))
})

test_that("infinite right-hand sides become bounds", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 2, lower = 0, upper = 5)
  opt_constraint(m, x <= c(Inf, 2), name = "a")
  opt_constraint(m, -Inf <= x[1], name = "b")
  opt_objective(m, -sum(x))
  sf <- opt_standard_form(m)
  expect_false(anyNA(c(sf$row_lower, sf$row_upper)))
  expect_equal(unname(opt_solve(m, solver = "highs")$x), c(5, 2))
  expect_error(opt_constraint(m, x[1] >= Inf), "never hold")
  z <- opt_variable(m, "z", type = "binary")
  expect_error(opt_indicator(m, z, x[1] <= Inf), "finite")
})
