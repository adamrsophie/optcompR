# Correctness of statuses, bounds, feasibility checks and result snapshots.
# Each case states the corrected behavior; most use tiny models whose
# answers are known exactly.
# A controlled backend exposes validator behavior independently of numerical
# solver behavior. Its deliberately wrong feasibility claim must not be trusted.
zero_backend_result <- function(model) {
  nm <- "test_zero_backend"
  opt_register_solver(nm, function(sf, cls, control) {
    list(termination = "optimal", primal_status = "feasible_point",
         dual_status = "feasible_point", x = rep(0, length(sf$c)),
         objective = 0, row_dual = rep(0, nrow(sf$A)),
         col_dual = rep(0, length(sf$c)), extra = list())
  }, classes = c("LP", "QP", "NLP"), version = "test double")
  on.exit(opt_unregister_solver(nm), add = TRUE)
  opt_solve(model, solver = nm)
}

capture_error <- function(expr) tryCatch(force(expr), error = identity)

# This helper preserves evidence if a future central result validator rejects
# the point before opt_check. A solver-boundary rejection must explain the
# numerical/feasibility problem, not an unrelated argument/runtime error.
expect_rejection_reason <- function(err) {
  expect_match(conditionMessage(err), "feasib|violat|finite|invalid.*point|invalid.*solution",
               ignore.case = TRUE)
}

make_soc_model <- function(integer = FALSE) {
  m <- opt_model()
  x <- opt_variable(m, "x")
  t <- opt_variable(m, "t")
  p <- opt_parameter(m, "p", 2)
  if (integer) {
    k <- opt_variable(m, "k", type = "integer", lower = 1, upper = 1)
    opt_constraint(m, x == p + k)
  } else opt_constraint(m, x == p)
  opt_soc_constraint(m, x, t)
  opt_objective(m, t)
  list(model = m, x = x, t = t, p = p)
}

expect_snapshot_derivative <- function(result, parameter) {
  dj <- capture_error(opt_solution_jacobian(result, parameter))
  expect_false(inherits(dj, "error"),
               info = if (inherits(dj, "error")) conditionMessage(dj) else NULL)
  if (!inherits(dj, "error")) {
    got <- dj$solution[c("x", "t"), , drop = FALSE]
    expect_equal(as.numeric(got), c(1, 1), tolerance = 1e-5)
  }
}

test_that("independent diagnostics retain upper-only and lower-only violations", {
  for (direction in c("lower", "upper")) {
    m <- opt_model()
    x <- opt_variable(m, "x")
    opt_objective(m, x^2)
    opt_constraint(m, if (direction == "lower") x >= 1 else x <= -1)
    r <- capture_error(zero_backend_result(m))
    if (inherits(r, "error")) {
      expect_rejection_reason(r)
    } else {
      # Force the diagnostic fixture to retain a point, even if central result
      # normalization now marks it infeasible. opt_check is an independent check.
      r$primal_status <- "feasible_point"
      d <- opt_check(r)
      expect_false(d$feasible, info = direction)
      expect_gte(d$max_row_violation, 1 - 1e-12)
      expect_equal(d$rows$violation, 1, info = direction)
    }
  }
})

test_that("solution derivatives reject infeasible stationary points", {
  for (where in c("row", "bound")) {
    m <- opt_model()
    x <- opt_variable(m, "x", lower = if (where == "bound") 1 else -Inf)
    p <- opt_parameter(m, "p", 0)
    if (where == "row") opt_constraint(m, x >= 1)
    opt_objective(m, (x - p)^2)
    r <- capture_error(zero_backend_result(m))
    if (inherits(r, "error")) {
      expect_rejection_reason(r)
    } else {
      r$primal_status <- "feasible_point"
      expect_error(opt_solution_jacobian(r, p),
                   "feasib|violat|KKT|differentiable", ignore.case = TRUE,
                   info = where)
    }
  }
})

test_that("NLopt does not certify an impossible exponential row", {
  skip_if_not_installed("nloptr")
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 0, upper = 2, start = 1)
  opt_constraint(m, exp(x) <= 0)
  opt_objective(m, x^2)
  r <- capture_error(opt_solve(m, solver = "nloptr",
    options = list(algorithm = "NLOPT_LN_COBYLA")))
  if (inherits(r, "error")) {
    expect_rejection_reason(r)
  } else {
    expect_false(opt_primal_status(r) == "feasible_point")
    expect_false(opt_status(r) %in% c("optimal", "locally_optimal", "almost_optimal"))
    expect_false(opt_has_solution(r))
  }
})

test_that("NLP branch-and-bound reports a valid bound and honest gap", {
  skip_if_not_installed("nloptr")
  m <- opt_model()
  x <- opt_variable(m, "x", type = "integer", lower = 0, upper = 1)
  opt_objective(m, (x - 0.6)^2)
  r <- opt_solve(m, solver = "nlpbb", mip_gap = 0.5)
  exact <- min((0:1 - 0.6)^2)
  incumbent <- opt_objective_value(r)
  bound <- opt_objective_bound(r)
  expect_true(is.finite(bound))
  expect_lte(bound, exact + 1e-6)
  expect_gte(incumbent, exact - 1e-6)
  if (opt_status(r) == "optimal")
    expect_equal(incumbent, exact, tolerance = 1e-6)
  if (incumbent > exact + 1e-6) {
    expect_identical(opt_status(r), "gap_limit")
    expect_true(is.finite(opt_gap(r)) && opt_gap(r) > 0)
  }
})

test_that("truncated SCS iterates are labeled by actual primal feasibility", {
  skip_if_not_installed("scs")
  m <- opt_model()
  x <- opt_variable(m, "x")
  opt_constraint(m, x >= 10)
  opt_objective(m, x)
  r <- opt_solve(m, solver = "scs", options = list(max_iters = 1L))
  # r$x is used directly: opt_value may correctly warn once the point status is
  # repaired. No exact iterate or runtime is asserted.
  genuinely_feasible <- length(r$x) == 1L && is.finite(r$x[1]) && r$x[1] >= 10 - 1e-5
  if (!genuinely_feasible) {
    expect_false(opt_primal_status(r) == "feasible_point")
    expect_false(opt_has_solution(r))
  } else {
    expect_true(is.finite(opt_objective_value(r)))
  }
})

test_that("infeasibility is invariant to positive penalty rescaling", {
  m <- opt_model()
  x <- opt_variable(m, "x")
  opt_constraint(m, x >= 1, name = "lower")
  opt_constraint(m, x <= 0, name = "upper")
  regular <- opt_explain_infeasibility(m, bounds = FALSE,
    weights = c(lower = 1, upper = 1), solver = "highs")
  tiny <- opt_explain_infeasibility(m, bounds = FALSE,
    weights = c(lower = 1e-9, upper = 1e-9), solver = "highs")
  expect_false(regular$feasible)
  expect_false(tiny$feasible)
  expect_identical(tiny$feasible, regular$feasible)
  expect_gte(max(tiny$violations$violation), 0.5)
})

test_that("conic helper results remain differentiable after restoration", {
  skip_if_not_installed("clarabel")
  b <- make_soc_model()
  direct <- opt_solve(b$model, solver = "clarabel")
  expect_snapshot_derivative(direct, b$p) # positive control
  weighted <- opt_weighted(b$model, list(b$t, 2 * b$t), weights = c(1, 1),
                           solver = "clarabel")
  expect_snapshot_derivative(weighted, b$p)
  b <- make_soc_model()
  ss <- opt_scenarios(b$model, b$p, c(2, 3), keep = TRUE, solver = "clarabel")
  expect_length(attr(ss, "results"), 2L)
  for (r in attr(ss, "results")) expect_snapshot_derivative(r, b$p)
})

test_that("fixed relaxation is reusable and returns a useful snapshot", {
  m <- opt_model()
  x <- opt_variable(m, "x", type = "integer", lower = 0, upper = 3)
  opt_objective(m, x)
  r <- opt_solve(m, solver = "highs")
  before <- opt_variables(m)
  fr1 <- opt_fixed_relaxation(r, solver = "highs")
  expect_identical(opt_variables(m), before)
  fr2 <- capture_error(opt_fixed_relaxation(r, solver = "highs"))
  expect_false(inherits(fr2, "error"),
               info = if (inherits(fr2, "error")) conditionMessage(fr2) else NULL)
  if (!inherits(fr2, "error"))
    expect_equal(opt_objective_value(fr2), opt_objective_value(fr1), tolerance = 1e-8)
  # Optional extension of the same contract: fixed MISOCP -> continuous SOCP.
  if (requireNamespace("scip", quietly = TRUE) &&
      requireNamespace("clarabel", quietly = TRUE)) {
    b <- make_soc_model(integer = TRUE)
    mi <- opt_solve(b$model, solver = "scip")
    fr <- opt_fixed_relaxation(mi, solver = "clarabel")
    expect_snapshot_derivative(fr, b$p)
  }
})

test_that("logical helpers reject non-Boolean affine literals", {
  for (kind in c("scaled", "sum", "bad_complement")) {
    m <- opt_model()
    a <- opt_variable(m, "a", type = "binary")
    b <- opt_variable(m, "b", type = "binary")
    bad <- switch(kind, scaled = 2 * a, sum = a + b, bad_complement = 1 - 2 * a)
    expect_error(opt_implies(m, bad, b), "literal|binary|Boolean|boolean",
                 info = kind)
  }
  m <- opt_model()
  a <- opt_variable(m, "a", type = "binary")
  b <- opt_variable(m, "b", type = "binary")
  expect_no_error(opt_implies(m, 1 - a, b))
})

test_that("rejected robust-constraint construction is atomic", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 2)
  before_vars <- opt_variables(m)
  before_cons <- opt_constraints(m)
  expect_error(opt_robust_constraint(m, x, a0 = 1, b = 1, delta = 1,
                                    uncertainty = "budget", gamma = 3), "gamma")
  expect_identical(opt_variables(m), before_vars)
  expect_identical(opt_constraints(m), before_cons)
})

test_that("SCIP polynomial lifting preserves solver objectives and bounds", {
  skip_if_not_installed("scip")
  for (offset in c(0, 5, -3)) {
    m <- opt_model()
    x <- opt_variable(m, "x", lower = -2, upper = 2)
    opt_objective(m, x^4 + offset)
    r <- opt_solve(m, solver = "scip")
    expect_equal(opt_objective_value(r), offset, tolerance = 1e-5)
    expect_true(is.finite(opt_objective_value(r, solver = TRUE)))
    expect_equal(opt_objective_value(r, solver = TRUE), offset, tolerance = 1e-5)
    expect_true(is.finite(opt_objective_bound(r)))
    expect_equal(opt_objective_bound(r), offset, tolerance = 1e-5)
  }
})

test_that("expression simplification and numeric matrix shapes are preserved", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 3)
  expect_equal(opt_degree(x - x + 1), rep(0L, 3))
  expect_equal(opt_degree(x^2 - x^2), rep(0L, 3))
  B <- matrix(1:3, nrow = 1)
  expr <- capture_error(x %*% B)
  expect_false(inherits(expr, "error"),
               info = if (inherits(expr, "error")) conditionMessage(expr) else NULL)
  if (!inherits(expr, "error"))
    expect_equal(opt_evaluate(expr, 1:3), (1:3) %*% B)
})

test_that("a backend declaring no global classes gives local optima on nonconvex models", {
  nm <- "test_local_backend"
  opt_register_solver(nm, function(sf, cls, control) {
    list(termination = "optimal", primal_status = "feasible_point",
         dual_status = "not_available", x = rep(0.5, length(sf$c)),
         objective = NA_real_)
  }, classes = c("NLP"), global = character(), version = "test double")
  on.exit(opt_unregister_solver(nm), add = TRUE)
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 0, upper = 1)
  opt_objective(m, sin(3 * x))
  expect_identical(opt_status(opt_solve(m, solver = nm)), "locally_optimal")
})

test_that("SCS dual infeasibility is not reported as unboundedness", {
  skip_if_not_installed("scs")
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 0)
  y <- opt_variable(m, "y")
  opt_constraint(m, y >= 1)
  opt_constraint(m, y <= 0)
  opt_objective(m, -x)
  st <- opt_status(opt_solve(m, solver = "scs"))
  expect_false(st == "unbounded")
})

test_that("a native gap option and mip_gap give the same status", {
  set.seed(80)
  n <- 80
  w <- runif(n, 1, 10)
  v <- w + runif(n, 0, 5)
  build <- function() {
    m <- opt_model("max")
    x <- opt_variable(m, "x", dim = n, type = "binary")
    opt_constraint(m, sum(w * x) <= sum(w) / 3)
    opt_objective(m, sum(v * x))
    m
  }
  a <- opt_solve(build(), solver = "highs", mip_gap = 0.5)
  b <- opt_solve(build(), solver = "highs", options = list(mip_rel_gap = 0.5))
  expect_identical(opt_status(a), opt_status(b))
  if (opt_gap(a) > 1e-9) expect_identical(opt_status(a), "gap_limit")
})

test_that("two-sided quadratic rows are classified by redundancy of a side", {
  conv <- function(f) {
    m <- opt_model()
    x <- opt_variable(m, "x")
    f(m, x)
    opt_objective(m, x)
    opt_classify(m)$convexity
  }
  expect_identical(conv(function(m, x) opt_constraint(m, x^2, lower = 0, upper = 1)), "convex")
  expect_identical(conv(function(m, x) opt_constraint(m, x^2 == 0)), "convex")
  expect_identical(conv(function(m, x) opt_constraint(m, x^2, lower = 1, upper = 4)), "nonconvex")
})

test_that("named points reject repeated names and wrong lengths", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 2)
  expect_error(opt_evaluate(sum(x), c("x[1]" = 5, "x[1]" = 100, "x[2]" = 7)),
               "one value per variable|repeat")
  expect_error(opt_evaluate(sum(x), c("x[1]" = 5, "x[1]" = 100)), "repeat")
  expect_equal(opt_evaluate(sum(x), c("x[2]" = 7, "x[1]" = 5)), 12)
})

test_that("the elastic relaxation minimizes total violation, not rows changed", {
  m <- opt_model()
  x <- opt_variable(m, "x")
  opt_constraint(m, x >= 1, name = "a")
  opt_constraint(m, x >= 1, name = "b")
  opt_constraint(m, 100 * x <= 0, name = "c")
  e <- opt_explain_infeasibility(m, bounds = FALSE, solver = "highs")
  expect_false(e$feasible)
  expect_setequal(e$violations$name, c("a", "b"))
  expect_equal(e$total, 2, tolerance = 1e-6)
})
