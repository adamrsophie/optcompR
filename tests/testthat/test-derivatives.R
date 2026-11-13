fd_grad <- function(f, x, h = 1e-6) {
  vapply(seq_along(x), function(i) {
    e <- replace(numeric(length(x)), i, h)
    (f(x + e) - f(x - e)) / (2 * h)
  }, 0)
}

test_that("QP gradients and Hessians are structural and exact", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 3)
  Q <- matrix(c(4, 1, 0, 1, 3, -1, 0, -1, 2), 3)
  cc <- c(1, -2, 0.5)
  opt_objective(m, 0.5 * (x %*% (Q %*% x)) + sum(cc * x))
  pt <- c(0.3, -1.2, 2)
  expect_equal(unname(opt_gradient(m, pt)), as.vector(Q %*% pt + cc))
  expect_equal(as.matrix(opt_hessian(m, pt)), Q, ignore_attr = TRUE)
  # asymmetric input matrix: gradient is 0.5 (Q + Q') x + c
  A <- matrix(c(1, 2, 0, 0, 1, 3, 1, 0, 2), 3)
  opt_objective(m, x %*% (A %*% x) + sum(cc * x))
  expect_equal(unname(opt_gradient(m, pt)), as.vector((A + t(A)) %*% pt + cc))
  expect_equal(as.matrix(opt_hessian(m, pt)), A + t(A), ignore_attr = TRUE)
})

test_that("nonlinear derivatives agree with central differences", {
  m <- opt_model()
  x <- opt_variable(m, "x")
  y <- opt_variable(m, "y")
  opt_objective(m, exp(x * y) + log(1 + x^2) * sin(y) + sqrt(x + 3) / (y^2 + 1) +
                  tanh(x - y) + atan(x)^2 + (x + 2)^1.5)
  f <- function(v) {
    x <- v[1]; y <- v[2]
    exp(x * y) + log(1 + x^2) * sin(y) + sqrt(x + 3) / (y^2 + 1) +
      tanh(x - y) + atan(x)^2 + (x + 2)^1.5
  }
  for (pt in list(c(0.2, 0.7), c(-0.5, 1.3))) {
    expect_equal(unname(opt_gradient(m, pt)), fd_grad(f, pt), tolerance = 1e-7)
    H <- as.matrix(opt_hessian(m, pt))
    Hfd <- rbind(fd_grad(function(v) opt_gradient(m, v)[1], pt),
                 fd_grad(function(v) opt_gradient(m, v)[2], pt))
    expect_equal(H, Hfd, tolerance = 1e-6, ignore_attr = TRUE)
  }
})

test_that("constraint Jacobians and Lagrangian Hessians", {
  m <- opt_model()
  x <- opt_variable(m, "x")
  y <- opt_variable(m, "y")
  opt_constraint(m, x * y <= 4, name = "prod")
  opt_constraint(m, c(exp(x), x + 2 * y) <= c(3, 1), name = "pair")
  opt_objective(m, x^2 + y^2)
  pt <- c(x = 0.5, y = 2)
  J <- opt_jacobian(m, pt)
  expect_equal(dim(J), c(3L, 2L))
  expect_equal(as.matrix(J), rbind(c(2, 0.5), c(exp(0.5), 0), c(1, 2)),
               ignore_attr = TRUE)
  expect_equal(rownames(J), c("prod", "pair[1]", "pair[2]"))
  H <- opt_hessian(m, pt, multipliers = c(2, 1, 0))
  # 2I - 2 * Hess(xy) - 1 * Hess(exp(x))
  expect_equal(as.matrix(H), rbind(c(2 - exp(0.5), -2), c(-2, 2)),
               ignore_attr = TRUE)
  expect_error(opt_gradient(m, c(1, 2, 3)), "one value per variable")
})

test_that("derivatives at a result use stored parameter values", {
  m <- opt_model()
  x <- opt_variable(m, "x")
  p <- opt_parameter(m, "p", value = 3)
  opt_objective(m, p * x^2)
  r <- opt_solve(m)
  expect_equal(unname(opt_gradient(m, c(x = 2))), 12)
  expect_equal(unname(opt_gradient(m, r)), 0, tolerance = 1e-8)
})

test_that("models are written to MPS and LP files", {
  m <- opt_model("max")
  x <- opt_variable(m, "x", dim = 2, lower = 0, upper = 3)
  y <- opt_variable(m, "y", type = "integer", lower = 0, upper = 2)
  opt_constraint(m, x[1] + 2 * x[2] + y <= 4)
  opt_objective(m, 3 * x[1] + x[2] + y)
  for (ext in c(".mps", ".lp")) {
    f <- tempfile(fileext = ext)
    opt_write_model(m, f)
    expect_true(file.exists(f))
    expect_gt(length(readLines(f)), 3)
  }
  expect_error(opt_write_model(m, tempfile(fileext = ".txt")), "mps")
  opt_objective(m, exp(x[1]))
  expect_error(opt_write_model(m, tempfile(fileext = ".mps")), "not supported")
})
