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

# nine parameters and one state: "auto" selects the adjoint
decay9 <- function(...) {
  modelcompR::ODEProblem(function(u, p, t) -sum(p) * u, 1, c(0, 1), rep(.1, 9), ...)
}
methods3 <- c("forward", "adjoint", "auto")

test_that("ODE gradients include parameter-dependent initial values", {
  skip_if_not_installed("modelcompR")
  skip_if_not_installed("symcompR")
  t0 <- modelcompR::t_nounits
  D <- modelcompR::D_nounits
  u <- modelcompR::Variable("tinit_u", t0, default = symcompR::Symbol("tinit_a") * 2)
  ks <- lapply(seq_len(8), function(j)
    modelcompR::Parameter(paste0("tinit_k", j), default = .5 / 8))
  a <- modelcompR::Parameter("tinit_a", default = 2)
  rate <- Reduce(`+`, ks)
  sys <- modelcompR::mtkcompile(modelcompR::System(
    list(D(u) ~ -rate * u + 0 * a), t0, name = "tinit"))
  prob <- modelcompR::ODEProblem(sys, list(), c(0, 1))
  v <- as.numeric(prob$p)
  j <- match("tinit_a", names(prob$p))
  m <- opt_model()
  theta <- opt_variable(m, "theta", dim = length(v))
  times <- c(.5, 1)
  # loss sum (2 a exp(-t/2))^2, derivative 8 a sum exp(-t) at a = 2
  for (method in methods3) {
    opt_objective(m, opt_ode_loss(prob, theta, times, matrix(0, 2, 1),
                                  sensitivity = method))
    expect_equal(unname(opt_gradient(m, v))[j], 16 * sum(exp(-times)),
                 tolerance = 1e-5, info = method)
  }
})

test_that("a state observed twice adds both contributions", {
  skip_if_not_installed("modelcompR")
  v <- rep(.1, 9)
  times <- c(.2, .5, 1)
  m <- opt_model()
  theta <- opt_variable(m, "theta", dim = 9)
  data <- cbind(exp(-.7 * times), exp(-1.2 * times))
  w <- matrix(rep(c(2, 3), each = 3), 3)
  u <- exp(-sum(v) * times)
  g1 <- sum(2 * w * (cbind(u, u) - data) * cbind(-times * u, -times * u))
  for (method in methods3) {
    opt_objective(m, opt_ode_loss(decay9(), theta, times, data, states = c(1, 1),
                                  weights = w, sensitivity = method))
    expect_equal(unname(opt_gradient(m, v))[1], g1, tolerance = 1e-6, info = method)
  }
})

test_that("ODE observations outside the interval are rejected", {
  skip_if_not_installed("modelcompR")
  m <- opt_model()
  theta <- opt_variable(m, "theta", dim = 9)
  for (method in methods3) {
    expect_error(opt_ode_loss(decay9(), theta, c(2, 3), matrix(0, 2, 1),
                              sensitivity = method), "interval")
  }
})

test_that("ODE losses do not depend on the order of observations", {
  skip_if_not_installed("modelcompR")
  v <- rep(.1, 9)
  times <- c(.2, .5, 1, .5)
  o <- c(3, 1, 4, 2)
  data <- matrix(exp(-.7 * times), ncol = 1)
  u <- exp(-sum(v) * times)
  m <- opt_model()
  theta <- opt_variable(m, "theta", dim = 9)
  for (method in methods3) {
    e <- opt_ode_loss(decay9(), theta, times[o], data[o, , drop = FALSE],
                      sensitivity = method)
    expect_equal(as.numeric(opt_evaluate(e, v)), sum((u - data)^2),
                 tolerance = 1e-8, info = method)
    opt_objective(m, e)
    expect_equal(unname(opt_gradient(m, v))[1],
                 sum(2 * (u - data) * (-times * u)), tolerance = 1e-6, info = method)
  }
})

test_that("every ODE gradient method supports callbacks", {
  skip_if_not_installed("modelcompR")
  v <- rep(.1, 9)
  cb <- modelcompR::DiscreteCallback(function(u, t, integrator) t == .5,
    function(integrator) integrator$u <- 2 * integrator$u)
  prob <- decay9(callback = cb, tstops = .5)
  times <- c(.2, .7, 1)
  # u(t) = exp(-s t) before 0.5 and 2 exp(-s t) after, with s = sum(p)
  s <- sum(v)
  u <- exp(-s * times) * ifelse(times > .5, 2, 1)
  m <- opt_model()
  theta <- opt_variable(m, "theta", dim = 9)
  for (method in methods3) {
    e <- opt_ode_loss(prob, theta, times, matrix(0, 3, 1), sensitivity = method)
    expect_equal(as.numeric(opt_evaluate(e, v)), sum(u^2), tolerance = 1e-7,
                 info = method)
    opt_objective(m, e)
    expect_equal(unname(opt_gradient(m, v))[1], sum(-2 * times * u^2),
                 tolerance = 1e-5, info = method)
  }
})
