# Construction paths that must stay linear in the size of the model, and
# the compact representations they produce.

test_that("c() stacks polynomial and nonlinear parts in order", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 3)
  p <- opt_parameter(m, "p", 2)
  e <- c(x[1], 2 * x[2]^2 + 1, exp(x[3]), 5, p * x[1])
  expect_identical(e$len, 5L)
  expect_equal(opt_evaluate(e, c(1, 2, 0)), c(1, 9, 1, 5, 2))
  xs <- lapply(1:2000, function(i) x[(i %% 3) + 1] + i)
  big <- do.call(c, xs)
  expect_equal(opt_evaluate(big, c(1, 2, 3)),
               c(1, 2, 3)[((1:2000) %% 3) + 1] + 1:2000)
  m2 <- opt_model()
  y <- opt_variable(m2, "y")
  expect_error(c(x[1], y), "different models")
})

test_that("quadratic rows are stored as one triplet table", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 3)
  opt_constraint(m, x^2 <= 1)
  opt_constraint(m, x[1] * x[2] + x[3]^2 <= 2)
  sf <- opt_standard_form(m)
  expect_setequal(unique(sf$row_quad$i), 1:4)
  H <- opt_row_hessian(sf, 4)
  expect_equal(as.matrix(H), matrix(c(0, 1, 0, 1, 0, 0, 0, 0, 2), 3))
  expect_null(opt_row_hessian(sf, 5))
  expect_lt(as.numeric(object.size(sf$row_quad)), 2000)
})

test_that("indicator rows are extracted without dense rows", {
  n <- 500
  m <- opt_model("max")
  x <- opt_variable(m, "x", dim = n, lower = 0, upper = 10)
  z <- opt_variable(m, "z", dim = n, type = "binary")
  opt_indicator(m, z, x + 0 * rev(x) <= 1, name = "switch")
  opt_objective(m, sum(x) - sum(z))
  sf <- opt_standard_form(m)
  expect_length(sf$indicators, n)
  expect_identical(sf$indicators[[3]]$j, 3L)
  expect_identical(sf$indicators[[3]]$rhs, 1)
  r <- opt_solve(m, solver = "highs")
  expect_equal(opt_objective_value(r), 10 * n)
})

test_that("kept HiGHS and OSQP workspaces give the answers of fresh solves", {
  build <- function(bv, qv) {
    m <- opt_model("max")
    x <- opt_variable(m, "x", dim = 3, lower = 0, upper = 4)
    b <- opt_parameter(m, "b", value = bv)
    q <- opt_parameter(m, "q", value = qv)
    opt_constraint(m, x[1] + 2 * x[2] + q * x[3] <= b)
    opt_objective(m, 3 * x[1] + 2 * x[2] + x[3])
    list(m = m, b = b, q = q)
  }
  k <- build(4, 1)
  r1 <- opt_solve(k$m, solver = "highs")
  for (v in list(c(6, 1), c(3, 2), c(5, 0.5))) {
    opt_set_parameter(k$m, k$b, v[1])
    opt_set_parameter(k$m, k$q, v[2])
    r <- opt_solve(k$m, solver = "highs")
    expect_true(r$extra$reused)
    fresh <- opt_solve(build(v[1], v[2])$m, solver = "highs")
    expect_equal(opt_objective_value(r), opt_objective_value(fresh),
                 tolerance = 1e-9)
  }
  # ranging is computed on request, from the stored basis
  s <- opt_sensitivity(r)
  expect_identical(nrow(s$costs), 3L)
  # a model read back from a file builds a new native solver
  f <- tempfile(fileext = ".rds")
  saveRDS(k$m, f)
  m2 <- readRDS(f)
  r2 <- opt_solve(m2, solver = "highs")
  expect_equal(opt_objective_value(r2), opt_objective_value(r), tolerance = 1e-9)
  skip_if_not_installed("osqp")
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 4, lower = -1, upper = 1)
  a <- opt_parameter(m, "a", dim = 4, value = c(0.5, 2, -3, 0))
  opt_objective(m, sum((x - a)^2))
  opt_solve(m, solver = "osqp")
  opt_set_parameter(m, a, c(-0.2, 0.3, 4, -5))
  r <- opt_solve(m, solver = "osqp")
  expect_true(r$extra$reused)
  expect_equal(unname(opt_value(r, x)), c(-0.2, 0.3, 1, -1), tolerance = 1e-5)
})
