bi_model <- function() {
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 0, upper = 4)
  y <- opt_variable(m, "y", lower = 0, upper = 4)
  opt_constraint(m, x + y >= 4)
  list(m = m, x = x, y = y)
}

test_that("weighted sums restore the model objective", {
  b <- bi_model()
  opt_objective(b$m, b$x)
  r <- opt_weighted(b$m, list(cost = 2 * b$x + b$y, risk = b$x + 3 * b$y),
                    c(1, 1))
  # minimize 3x + 4y on x + y >= 4: x = 4, y = 0
  expect_equal(unname(r$objective_values), c(8, 4))
  expect_identical(names(r$objective_values), c("cost", "risk"))
  expect_equal(optcompR:::eval_expr(b$m$objective, c(5, 0)), 5)
  expect_error(opt_weighted(b$m, list(b$x, b$y), c(1, -1)), "nonnegative")
})

test_that("lexicographic stages respect priorities", {
  b <- bi_model()
  st <- opt_lexicographic(b$m, list(cost = b$x + b$y, risk = b$x))
  expect_length(st, 2)
  expect_equal(opt_value(st[[2]], c(b$x, b$y)), c(0, 4), tolerance = 1e-6)
  expect_equal(nrow(opt_constraints(b$m)), 1)
})

test_that("epsilon constraints trace the front and clean up", {
  b <- bi_model()
  front <- opt_epsilon_constraint(b$m, 2 * b$x + b$y, b$x + 3 * b$y,
                                  epsilon = c(4, 8, 12, 2))
  expect_equal(front$status, c("optimal", "optimal", "optimal", "infeasible"))
  # on the edge x + y = 4: primary = 4 + x, secondary = 12 - 2x
  expect_equal(front$primary[1:3], c(8, 6, 4), tolerance = 1e-7)
  expect_equal(nrow(opt_constraints(b$m)), 1)
})
