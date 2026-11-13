test_that("ODE parameters are recovered with exact sensitivities", {
  skip_if_not_installed("modelcompR")
  skip_if_not_installed("nloptr")
  f <- function(u, p, t) c(-p[1] * u[1], p[1] * u[1] - p[2] * u[2])
  prob <- modelcompR::ODEProblem(f, c(1, 0), c(0, 3), c(0.5, 0.5))
  times <- c(0.5, 1, 1.5, 2, 3)
  truth <- c(0.7, 0.3)
  # analytic solution of the two-compartment chain
  u1 <- exp(-truth[1] * times)
  u2 <- truth[1] / (truth[2] - truth[1]) * (exp(-truth[1] * times) -
                                               exp(-truth[2] * times))
  y <- cbind(u1, u2)
  m <- opt_model()
  th <- opt_variable(m, "theta", dim = 2, lower = 0.01, upper = 3,
                     start = c(0.4, 0.6))
  loss <- opt_ode_loss(prob, th, times, y)
  opt_objective(m, loss)
  # gradient matches central differences of the loss
  g <- opt_gradient(m, c(0.6, 0.4))
  val <- function(v) optcompR:::eval_expr(loss, c(v))
  h <- 1e-5
  fd <- vapply(1:2, function(j) {
    e <- replace(c(0, 0), j, h)
    (val(c(0.6, 0.4) + e) - val(c(0.6, 0.4) - e)) / (2 * h)
  }, 0)
  expect_equal(unname(g), fd, tolerance = 1e-5)
  r <- opt_solve(m)
  expect_equal(opt_value(r, th), truth, tolerance = 1e-5)
  expect_lt(opt_objective_value(r), 1e-10)
  expect_error(opt_ode_loss(prob, th, times, y[-1, ]), "one row per")
})
