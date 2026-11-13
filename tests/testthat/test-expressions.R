test_that("linear expressions carry exact coefficients", {
  m <- opt_model()
  x <- opt_variable(m, "x")
  y <- opt_variable(m, "y")
  e <- 3 * x - 2 * y + 5 - x
  z <- c(2, 7)
  expect_equal(optcompR:::eval_expr(e, z), 3 * 2 - 2 * 7 + 5 - 2)
  e <- optcompR:::compact_poly(e)
  expect_equal(sort(e$lv), c(-2, 2))
  expect_equal(e$const, 5)
})

test_that("quadratic products and powers expand correctly", {
  m <- opt_model()
  x <- opt_variable(m, "x")
  y <- opt_variable(m, "y")
  e <- (x - 1)^2 + (x + y) * (2 * y - 3)
  f <- function(x, y) (x - 1)^2 + (x + y) * (2 * y - 3)
  for (pt in list(c(0, 0), c(1.5, -2), c(-3, 4))) {
    expect_equal(optcompR:::eval_expr(e, pt), f(pt[1], pt[2]))
  }
  expect_identical(optcompR:::expr_kind_label(e), "quadratic")
})

test_that("cancelling terms reduce the degree", {
  m <- opt_model()
  x <- opt_variable(m, "x")
  e <- x * x - x^2 + 2 * x
  expect_identical(optcompR:::expr_kind_label(e), "affine")
  expect_identical(optcompR:::expr_kind_label(x - x + 1), "constant")
})

test_that("vector expressions recycle only length one", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 3)
  expect_length(x + 1, 3)
  expect_length(x + c(1, 2, 3), 3)
  expect_error(x + c(1, 2), "Non-conformable")
  y <- opt_variable(m, "y", dim = 2)
  expect_error(x + y, "Non-conformable")
})

test_that("sum, indexing and concatenation", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 4)
  z <- c(1, 2, 3, 4)
  expect_equal(optcompR:::eval_expr(sum(x), z), 10)
  expect_equal(optcompR:::eval_expr(x[c(4, 1)], z), c(4, 1))
  expect_equal(optcompR:::eval_expr(x[-1], z), c(2, 3, 4))
  expect_equal(optcompR:::eval_expr(x[c(TRUE, FALSE)], z), c(1, 3))
  expect_equal(optcompR:::eval_expr(c(x[2], 2 * x[3], 7), z), c(2, 6, 7))
  expect_equal(optcompR:::eval_expr(cumsum(x), z), cumsum(z))
  expect_equal(optcompR:::eval_expr(sum(x, x[1]), z), 11)
})

test_that("matrix variables, names and transposes", {
  m <- opt_model()
  X <- opt_variable(m, "X", dim = c(2, 3))
  expect_equal(dim(X), c(2L, 3L))
  z <- 1:6
  expect_equal(optcompR:::eval_expr(X[2, ], z), c(2, 4, 6))
  expect_equal(optcompR:::eval_expr(t(X), z), as.vector(t(matrix(1:6, 2))))
  expect_equal(dim(t(X)), c(3L, 2L))
  W <- opt_variable(m, "W", dimnames = list(c("a", "b"), c("u", "v")))
  zz <- c(z, 10, 20, 30, 40)
  expect_equal(optcompR:::eval_expr(W["b", "u"], zz), 20)
  v <- opt_variable(m, "v", dimnames = c("p", "q"))
  expect_equal(optcompR:::eval_expr(v["q"], c(zz, 5, 6)), 6)
})

test_that("matrix products with numeric and sparse matrices", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 3)
  A <- matrix(c(1, 2, 0, -1, 3, 4), 2)
  z <- c(1, -2, 5)
  expect_equal(optcompR:::eval_expr(A %*% x, z), as.vector(A %*% z))
  S <- Matrix::Matrix(A, sparse = TRUE)
  expect_equal(optcompR:::eval_expr(S %*% x, z), as.vector(A %*% z))
  expect_equal(optcompR:::eval_expr(c(1, 2, 3) %*% x, z), sum(c(1, 2, 3) * z))
  expect_equal(optcompR:::eval_expr(x %*% x, z), sum(z^2))
  X <- opt_variable(m, "X", dim = c(2, 2))
  zz <- c(z, 1, 2, 3, 4)
  B <- matrix(c(2, 0, 1, 1), 2)
  expect_equal(optcompR:::eval_expr(X %*% B, zz),
               as.vector(matrix(1:4, 2) %*% B))
  expect_equal(optcompR:::eval_expr(X %*% X, zz),
               as.vector(matrix(1:4, 2) %*% matrix(1:4, 2)))
})

test_that("comparisons build relations and never logical values", {
  m <- opt_model()
  x <- opt_variable(m, "x")
  r <- x + 1 <= 3
  expect_s3_class(r, "opt_relation")
  expect_s3_class(3 >= x, "opt_relation")
  expect_error(x < 3, "Strict")
  expect_error(x > 3, "Strict")
  expect_error(x != 3, "closed")
  expect_error(if (x <= 3) 1, "argument")
  expect_error(x & x, "Logical")
})

test_that("expressions from different models cannot be combined", {
  m1 <- opt_model()
  m2 <- opt_model()
  x <- opt_variable(m1, "x")
  y <- opt_variable(m2, "y")
  expect_error(x + y, "different models")
  expect_error(opt_constraint(m2, x <= 1), "different model")
  expect_error(opt_objective(m2, x), "different model")
})

test_that("unsupported reductions explain the alternative", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 2)
  expect_error(max(x), "epigraph")
  expect_error(x / 0, "zero")
})

test_that("printing shows readable algebra", {
  m <- opt_model()
  x <- opt_variable(m, "x")
  y <- opt_variable(m, "y")
  expect_output(print(2 * x - y + 3), "2\\*x - y \\+ 3")
  expect_output(print(x + y <= 4), "x \\+ y <= 4")
  expect_output(print(exp(x)), "exp\\(x\\)")
})
