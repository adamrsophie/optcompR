test_that("parameter exponents stay power nodes until values are known", {
  m <- opt_model()
  x <- opt_variable(m, "x", lower = -5, upper = 5)
  p <- opt_parameter(m, "p", value = 2)
  cc <- opt_parameter(m, "c", value = 3)
  e <- x^p
  expect_equal(e$op, "powp")
  opt_objective(m, e + cc * x)
  expect_equal(opt_classify(m)$class, "QP")
  r <- opt_solve(m, solver = "nloptr")
  expect_equal(unname(opt_value(r, x)), -1.5, tolerance = 1e-5)
  expect_equal(unname(opt_solution_jacobian(r, cc)$solution[1, 1]), -0.5,
               tolerance = 1e-5)
  expect_error(opt_solution_jacobian(r, p), "not differentiable in p")
  opt_set_parameter(m, p, 4)
  expect_equal(opt_classify(m)$class, "NLP")
})

test_that("derivatives of a^p in the exponent are exact for a positive base", {
  m <- opt_model()
  y <- opt_variable(m, "y", lower = 0.5, upper = 3)
  q <- opt_parameter(m, "q", value = 1.5)
  opt_objective(m, (y^q - 2)^2)
  r <- opt_solve(m, solver = "nloptr")
  y0 <- 2^(1 / 1.5)
  expect_equal(unname(opt_value(r, y)), y0, tolerance = 1e-5)
  expect_equal(unname(opt_solution_jacobian(r, q)$solution[1, 1]),
               -log(y0) / 1.5 * y0, tolerance = 1e-4)
})

test_that("parameter powers round trip through symcompR", {
  skip_if_not_installed("symcompR")
  m <- opt_model()
  x <- opt_variable(m, "x")
  p <- opt_parameter(m, "p", value = 3)
  s <- opt_to_sym(x^p)
  back <- opt_from_sym(m, s)
  expect_equal(back$op, "powp")
  expect_equal(format(back), format(x^p))
})

test_that("opt_copy() makes an independent model", {
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 0)
  opt_constraint(m, x <= 10, name = "cap")
  y <- opt_variable(m, "y", lower = 0, upper = 1)
  opt_sos1(m, c(x, y))
  opt_objective(m, x)
  m2 <- opt_copy(m, name = "copy")
  x2 <- opt_get(m2, "x")
  expect_true(identical(x2$model, m2))
  expect_true(identical(m2$cons[[1]]$expr$model, m2))
  opt_set_bounds(m2, x2, lower = 2)
  opt_delete_constraint(m2, "cap")
  expect_equal(unname(opt_objective_value(opt_solve(m2))), 2)
  expect_equal(nrow(opt_constraints(m)), 1L)
  expect_equal(unname(opt_objective_value(opt_solve(m))), 0)
  expect_equal(m2$name, "copy")
})

test_that("special constraints can be listed and deleted", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 3, lower = 0, upper = 1)
  z <- opt_variable(m, "z", type = "binary")
  s1 <- opt_sos1(m, x, name = "pick")
  opt_indicator(m, z, sum(x) <= 0.5, name = "ind")
  sp <- opt_special_constraints(m)
  expect_equal(sp$name, c("pick", "ind"))
  expect_equal(sp$type, c("sos1", "indicator"))
  expect_equal(sp$size, c(3L, 1L))
  opt_delete_constraint(m, s1)
  expect_equal(opt_special_constraints(m)$name, "ind")
  opt_delete_constraint(m, opt_special_constraint(m, "indicator", "ind"))
  expect_equal(nrow(opt_special_constraints(m)), 0L)
  expect_null(m$special)
  expect_error(opt_special_constraint(m, "sos", "pick"), "No sos constraint")
  m2 <- opt_model()
  expect_error(opt_delete_constraint(m2, s1), "different model")
})

test_that("expression introspection reports terms, degrees and symbols", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 2)
  p <- opt_parameter(m, "p", value = 2)
  e <- c(3 * x[1] + 2, x[1] * x[2] - p * x[2])
  tr <- opt_terms(e)
  expect_equal(tr$constant, c(2, 0))
  expect_equal(nrow(tr$linear), 1L)
  expect_equal(tr$quadratic$coef, c(1, -1))
  expect_false(tr$nonlinear)
  expect_true(opt_terms(exp(x[1]))$nonlinear)
  expect_equal(opt_degree(c(e, exp(x[1]), as_expr(1) + 0 * x[1])),
               c(1L, 2L, NA, 0L))
  expect_setequal(opt_symbols(e[2]), c("x[1]", "x[2]", "p"))
  expect_equal(opt_evaluate(x[1]^2 + p * x[2], c(2, 1)), 6)
  expect_equal(opt_evaluate(x, c("x[2]" = 5, "x[1]" = 4)), c(4, 5))
})

test_that("a stop at the user MIP gap has its own status", {
  set.seed(3)
  n <- 40
  w <- sample(20:60, n, TRUE)
  v <- w + sample(0:9, n, TRUE)
  m <- opt_model("max")
  x <- opt_variable(m, "x", dim = n, type = "binary")
  for (k in 1:3) opt_constraint(m, sum(sample(1:5, n, TRUE) * x) <= 3 * n)
  opt_constraint(m, sum(w * x) <= sum(w) / 2)
  opt_objective(m, sum(v * x))
  r <- opt_solve(m, solver = "highs", mip_gap = 0.2)
  if (opt_gap(r) > 1e-9) {
    expect_equal(opt_status(r), "gap_limit")
    expect_equal(opt_primal_status(r), "feasible_point")
  } else {
    expect_equal(opt_status(r), "optimal")
  }
  expect_equal(opt_status(opt_solve(m, solver = "highs")), "optimal")
})

test_that("SCIP dual bounds are recovered from the reported gap", {
  expect_equal(scip_dual_bound(110, 0.1, "min"), 100)
  expect_equal(scip_dual_bound(-100, 0.1, "min"), -110)
  expect_equal(scip_dual_bound(100, 0.1, "max"), 110)
  expect_equal(scip_dual_bound(-110, 0.1, "max"), -100)
  expect_equal(scip_dual_bound(5, 0, "min"), 5)
  expect_true(is.na(scip_dual_bound(5, 1e20, "min")))
  skip_if_not_installed("scip")
  m <- opt_model("max")
  x <- opt_variable(m, "x", dim = 2, type = "integer", lower = 0, upper = 10)
  opt_constraint(m, 2 * x[1] + 3 * x[2] <= 12.5)
  opt_objective(m, x[1] + x[2])
  r <- opt_solve(m, solver = "scip")
  expect_equal(unname(opt_objective_bound(r)), unname(opt_objective_value(r)))
})

test_that("user backends can be registered, used and removed", {
  calls <- 0L
  lp <- function(sf, cls, control) {
    calls <<- calls + 1L
    s <- if (sf$sense == "max") -1 else 1
    r <- highs::highs_solve(L = s * sf$c, lower = sf$lower, upper = sf$upper,
                            A = sf$A, lhs = sf$row_lower, rhs = sf$row_upper)
    list(termination = "optimal", primal_status = "feasible_point",
         dual_status = "not_available", x = r$primal_solution,
         objective = s * r$objective_value + sf$offset)
  }
  opt_register_solver("test_lp", lp, classes = "LP", package = "highs",
                      priority = 1L)
  on.exit(opt_unregister_solver("test_lp"))
  expect_true("test_lp" %in% opt_solvers()$solver)
  m <- opt_model("max")
  x <- opt_variable(m, "x", lower = 0, upper = 3)
  opt_constraint(m, x <= 2)
  opt_objective(m, 2 * x + 1)
  r <- opt_solve(m)
  expect_equal(r$solver, "test_lp")
  expect_equal(unname(opt_objective_value(r)), 5)
  expect_true(opt_check(r)$feasible)
  expect_true(opt_check(r)$kkt$pass)
  bad <- function(sf, cls, control) list(termination = "done")
  opt_register_solver("test_bad", bad, classes = "LP")
  on.exit(opt_unregister_solver("test_bad"), add = TRUE)
  expect_error(opt_solve(m, solver = "test_bad"), "unknown termination")
  expect_error(opt_register_solver("x", lp, classes = "LPX"), "subset")
  expect_error(opt_register_solver("x", lp, classes = "LP",
                                   features = list(fast = TRUE)),
               "Unknown feature")
  expect_true(opt_unregister_solver("test_bad"))
  expect_false(opt_unregister_solver("test_bad"))
})

test_that("SCS agrees with Clarabel on cone programs and their duals", {
  skip_if_not_installed("scs")
  skip_if_not_installed("clarabel")
  tight <- list(eps_abs = 1e-9, eps_rel = 1e-9, max_iters = 100000L)
  both <- function(build) {
    a <- opt_solve(build(), solver = "clarabel")
    b <- opt_solve(build(), solver = "scs", options = tight)
    expect_equal(opt_status(b), "optimal")
    expect_equal(opt_objective_value(b), opt_objective_value(a),
                 tolerance = 1e-5)
    expect_equal(unname(b$row_dual), unname(a$row_dual), tolerance = 1e-4)
    expect_equal(lapply(b$extra$cone_dual, unname),
                 lapply(a$extra$cone_dual, unname), tolerance = 1e-4)
  }
  both(function() {
    m <- opt_model()
    x <- opt_variable(m, "x", dim = 2)
    t <- opt_variable(m, "t")
    opt_constraint(m, sum(x) == 2)
    opt_soc_constraint(m, x, t)
    opt_objective(m, t)
    m
  })
  both(function() {
    m <- opt_model("max")
    x <- opt_variable(m, "x")
    y <- opt_variable(m, "y")
    opt_constraint(m, y <= 1)
    opt_exp_cone(m, x, y, 2)
    opt_objective(m, x)
    m
  })
  both(function() {
    m <- opt_model()
    X <- opt_sym_variable(m, "X", 3)
    opt_psd_constraint(m, X)
    C <- matrix(c(2, 1, 0, 1, 2, 1, 0, 1, 2), 3)
    opt_constraint(m, X[1, 1] + X[2, 2] + X[3, 3] == 1)
    opt_objective(m, sum(C * X))
    m
  })
  expect_equal(psd_lower_order(3), c(1L, 2L, 4L, 3L, 5L, 6L))
})

test_that("opt_check() verifies conic KKT conditions with cone multipliers", {
  skip_if_not_installed("clarabel")
  m <- opt_model("max")
  x <- opt_variable(m, "x", lower = 0)
  y <- opt_variable(m, "y", lower = 0)
  z <- opt_variable(m, "z")
  opt_constraint(m, x + y <= 2)
  opt_pow_cone(m, x, y, z, 0.3, name = "pc")
  opt_objective(m, z)
  r <- opt_solve(m, solver = "clarabel")
  k <- opt_check(r)$kkt
  expect_true(k$pass)
  expect_equal(k$cones$cone, "pc")
  bad <- r
  bad$extra$cone_dual$pc <- -bad$extra$cone_dual$pc
  expect_false(opt_check(bad)$kkt$pass)
  m2 <- opt_model()
  X <- opt_sym_variable(m2, "X", 2)
  opt_psd_constraint(m2, X, name = "P")
  opt_constraint(m2, X[1, 1] + X[2, 2] == 1)
  opt_objective(m2, 2 * X[1, 1] + X[2, 2] + 2 * X[1, 2])
  r2 <- opt_solve(m2, solver = "clarabel")
  expect_true(opt_check(r2)$kkt$pass)
})

test_that("solutions of conic programs are differentiated through cone boundaries", {
  skip_if_not_installed("clarabel")
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 2)
  t <- opt_variable(m, "t")
  p <- opt_parameter(m, "p", value = 2)
  opt_constraint(m, sum(x) == p)
  opt_soc_constraint(m, x, t)
  opt_objective(m, t)
  d <- opt_solution_jacobian(opt_solve(m, solver = "clarabel"), p)
  expect_equal(unname(d$solution[, 1]), c(0.5, 0.5, 1 / sqrt(2)),
               tolerance = 1e-5)
  expect_equal(unname(d$objective), 1 / sqrt(2), tolerance = 1e-5)

  m <- opt_model("max")
  x <- opt_variable(m, "x")
  y <- opt_variable(m, "y")
  q <- opt_parameter(m, "q", value = 2)
  opt_constraint(m, y <= 1)
  opt_exp_cone(m, x, y, q)
  opt_objective(m, x)
  d <- opt_solution_jacobian(opt_solve(m, solver = "clarabel"), q)
  expect_equal(unname(d$solution["x", 1]), exp(-1), tolerance = 1e-5)

  m <- opt_model("max")
  x <- opt_variable(m, "x", lower = 0)
  y <- opt_variable(m, "y", lower = 0)
  z <- opt_variable(m, "z")
  b <- opt_parameter(m, "b", value = 2)
  opt_constraint(m, x + y <= b)
  opt_pow_cone(m, x, y, z, 0.3)
  opt_objective(m, z)
  d <- opt_solution_jacobian(opt_solve(m, solver = "clarabel"), b)
  expect_equal(unname(d$solution[, 1]), c(0.3, 0.7, 0.3^0.3 * 0.7^0.7),
               tolerance = 1e-4)

  m <- opt_model()
  X <- opt_sym_variable(m, "X", 2)
  w <- opt_parameter(m, "w", value = 1)
  opt_psd_constraint(m, X, name = "P")
  opt_constraint(m, X[1, 1] + X[2, 2] == w)
  opt_objective(m, 2 * X[1, 1] + X[2, 2] + 2 * X[1, 2])
  r <- opt_solve(m, solver = "clarabel")
  d <- opt_solution_jacobian(r, w)
  # X*(w) = w * v v' for the eigenvector v of the smallest eigenvalue
  e <- eigen(matrix(c(2, 1, 1, 1), 2), symmetric = TRUE)
  v <- e$vectors[, 2]
  V <- tcrossprod(v)
  expect_equal(unname(d$solution[, 1]), V[upper.tri(V, diag = TRUE)],
               tolerance = 1e-5)
  expect_equal(unname(d$objective), e$values[2], tolerance = 1e-5)
  r$extra$cone_dual <- NULL
  expect_error(opt_solution_jacobian(r, w), "cone multipliers")
})

test_that("active PSD derivatives match finite differences with parameters in the cone", {
  skip_if_not_installed("clarabel")
  build <- function(wv, pv) {
    m <- opt_model()
    X <- opt_sym_variable(m, "X", 3)
    w <- opt_parameter(m, "w", value = wv)
    p <- opt_parameter(m, "p", value = pv)
    C <- matrix(c(2, 1, 0, 1, 2, 1, 0, 1, 3), 3)
    opt_psd_constraint(m, X - p * diag(3))
    opt_constraint(m, X[1, 1] + X[2, 2] + X[3, 3] == w)
    opt_objective(m, sum(C * X) + X[1, 2])
    list(m = m, w = w, p = p)
  }
  b <- build(1, 0.1)
  r <- opt_solve(b$m, solver = "clarabel")
  d <- opt_solution_jacobian(r, b$p)
  h <- 1e-4
  r1 <- opt_solve(build(1, 0.1 + h)$m, solver = "clarabel")
  r0 <- opt_solve(build(1, 0.1 - h)$m, solver = "clarabel")
  expect_equal(unname(d$solution[, 1]), unname((r1$x - r0$x) / (2 * h)),
               tolerance = 1e-4)
  expect_equal(unname(d$objective),
               (opt_objective_value(r1) - opt_objective_value(r0)) / (2 * h),
               tolerance = 1e-4)
})

test_that("lowering is reused across parameter updates", {
  m <- opt_model("max")
  x <- opt_variable(m, "x", dim = 2, lower = 0)
  p <- opt_parameter(m, "p", value = 1)
  opt_constraint(m, x[1] + x[2] <= 4, name = "fixed")
  opt_constraint(m, x[1] <= 2 * p, name = "moving")
  opt_objective(m, x[1] + 2 * x[2])
  sf1 <- opt_standard_form(m)
  lc1 <- m$cache$lowered
  opt_set_parameter(m, p, 1.5)
  sf2 <- opt_standard_form(m)
  lc2 <- m$cache$lowered
  expect_identical(lc2$cons[[1]], lc1$cons[[1]])
  expect_false(identical(lc2$cons[[2]]$pv, lc1$cons[[2]]$pv))
  expect_equal(sf2$row_upper, c(4, 3))
  expect_equal(unname(opt_value(opt_solve(m), x)), c(0, 4))
  opt_constraint(m, x[2] <= 3, name = "new")
  expect_equal(opt_standard_form(m)$row_upper, c(4, 3, 3))
  expect_equal(m$cache$lowered$structure, m$structure_version)
})

test_that("scalar fast paths keep expression semantics", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 3)
  expect_equal(format(2 * x[2]), format(x[2] * 2))
  expect_equal(opt_terms(3 * x[3])$linear$coef, 3)
  expect_equal(opt_evaluate(x[2] * 2 + x[3], c(1, 2, 3)), 7)
  y <- x[2:3]
  expect_equal(opt_symbols(y[2]), "x[3]")
})

test_that("HiGHS restarts LPs from the previous basis after parameter updates", {
  set.seed(7)
  n <- 60
  A <- matrix(runif(40 * n), 40)
  m <- opt_model("max")
  x <- opt_variable(m, "x", dim = n, lower = 0, upper = 1)
  b <- opt_parameter(m, "b", dim = 40, value = rep(10, 40))
  opt_constraint(m, A %*% x <= b)
  opt_objective(m, sum(runif(n) * x))
  r1 <- opt_solve(m, solver = "highs")
  expect_false(r1$warm_start)
  opt_set_parameter(m, b, rep(10.5, 40))
  r2 <- opt_solve(m, solver = "highs")
  expect_true(r2$warm_start)
  cold <- opt_solve(opt_copy(m), solver = "highs")
  expect_false(cold$warm_start)
  expect_equal(opt_objective_value(r2), opt_objective_value(cold),
               tolerance = 1e-8)
  expect_lte(r2$iterations[["simplex"]], cold$iterations[["simplex"]])
  r3 <- opt_solve(m, solver = "highs", warm_start = FALSE)
  expect_false(r3$warm_start)
})

test_that("sparsity patterns contain every Jacobian and Hessian nonzero", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 4)
  p <- opt_parameter(m, "p", value = 1.5)
  fx <- opt_function(function(v) sum(v^2) * v[1], x[2:4], name = "f",
                    gradient = function(v) c(sum(v^2) + 2 * v[1]^2,
                                             2 * v[1] * v[-1]),
                    hessian = function(v) {
                      H <- 2 * v[1] * diag(length(v))
                      H[1, ] <- H[1, ] + 2 * v
                      H[, 1] <- H[, 1] + 2 * v
                      H[1, 1] <- 6 * v[1]
                      H
                    })
  es <- list(c(exp(x[1]) * x[2], x[3] / (1 + x[4]^2), sin(x[2]) + x[1]^3,
               x[4]^p, sum(x^2), fx),
             cumsum(log1p(x^2)), 3 * x + 1)
  ns <- n_symbols(m)
  set.seed(1)
  z <- c(runif(4, 0.5, 2), 1.5)
  for (e in es) {
    J <- as.matrix(jac_expr(e, z, ns)$J) != 0
    PJ <- as.matrix(jac_pattern(e, ns)) > 0
    expect_true(all(PJ[J]))
    H <- as.matrix(hess_expr(e, z, ns, runif(e$len) + 0.5)) != 0
    PH <- as.matrix(hess_pattern(e, ns)) > 0
    expect_true(all(PH[H]))
  }
  # separable sums stay diagonal
  PH <- as.matrix(hess_pattern(sum(exp(x)), ns))
  expect_equal(sum(PH > 0), 4L)
})

test_that("IPOPT solves a large sparse NLP with sparse structures", {
  skip_if_not(opt_ipopt_available(), "optcompR was built without IPOPT")
  n <- 400
  m <- opt_model()
  x <- opt_variable(m, "x", dim = n, lower = -5, upper = 5)
  opt_constraint(m, x[-1] - x[-n]^2 == 0)
  opt_objective(m, sum((x - 1)^2) + sum(exp(x / 10)))
  r <- opt_solve(m, solver = "ipopt")
  expect_true(opt_status(r) %in% c("optimal", "locally_optimal"))
  expect_true(opt_check(r)$kkt$pass)
})

test_that("composition rules prove convexity of nonlinear models", {
  cls <- function(build) {
    m <- opt_model()
    x <- opt_variable(m, "x", dim = 2, lower = 0.5, upper = 3)
    build(m, x)
    opt_classify(m)$convexity
  }
  expect_equal(cls(function(m, x) opt_objective(m, sum(exp(x)) - log(x[1]))),
               "convex")
  expect_equal(cls(function(m, x) {
    opt_objective(m, x[1] + x[2])
    opt_constraint(m, exp(x[1] + 2 * x[2]) <= 10)
    opt_constraint(m, sqrt(x[1]) + log(x[2]) >= 0.2)
    opt_constraint(m, 1 / x[1] + x[2]^2 <= 6)
  }), "convex")
  # concave objective for maximization
  expect_equal(cls(function(m, x) {
    opt_objective(m, log(x[1]) + sqrt(x[2]), sense = "max")
  }), "convex")
  # not provable: nonconvex or unknown, never "convex"
  expect_equal(cls(function(m, x) opt_objective(m, sin(x[1] * x[2]))),
               "unknown")
  expect_equal(cls(function(m, x) opt_objective(m, -exp(x[1]))), "unknown")
  expect_equal(cls(function(m, x) {
    opt_objective(m, x[1])
    opt_constraint(m, exp(x[1]) == 2)
  }), "unknown")
  # the bound domain matters: sin is concave on [0, pi]
  expect_equal(cls(function(m, x) opt_objective(m, -sin(x[1]), sense = "min")),
               "convex")
  m <- opt_model()
  y <- opt_variable(m, "y", lower = -4, upper = 4)
  opt_objective(m, -sin(y))
  expect_equal(opt_classify(m)$convexity, "unknown")
})

test_that("curvature intervals bound expression values", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 2, lower = c(-1, 0.5), upper = c(2, 3))
  es <- list(exp(x[1]) * 2 - x[2]^2, log(x[2]) + abs(x[1]), x[1]^3 / x[2],
             sqrt(x[2]) * x[1], cosh(x[1]) + tanh(x[2]))
  set.seed(3)
  pts <- cbind(runif(200, -1, 2), runif(200, 0.5, 3))
  ns <- n_symbols(m)
  slo <- c(-1, 0.5)
  shi <- c(2, 3)
  for (e in es) {
    cv <- curvature(e, slo, shi)
    v <- apply(pts, 1, function(z) eval_expr(e, z))
    expect_true(all(v >= cv$lo - 1e-12 & v <= cv$hi + 1e-12))
  }
})

test_that("NLP branch and bound certifies MINLPs proven convex by composition", {
  skip_if_not_installed("nloptr")
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 0, upper = 5)
  n <- opt_variable(m, "n", type = "integer", lower = 0, upper = 5)
  opt_objective(m, exp(x / 2) + (n - 2.6)^2 - x)
  opt_constraint(m, x + n >= 3)
  r <- opt_solve(m)
  expect_equal(r$classification$convexity, "convex")
  expect_equal(opt_status(r), "optimal")
  expect_equal(unname(opt_value(r, n)), 3)
})

test_that("infeasibility explanations name the rows and bounds that must move", {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 2, lower = 0, upper = c(1, 5))
  opt_constraint(m, x[1] + x[2] <= 1, name = "capacity")
  opt_constraint(m, x[1] >= 2, name = "demand_1")
  opt_constraint(m, x[2] >= 0.5, name = "demand_2")
  e <- opt_explain_infeasibility(m)
  expect_false(e$feasible)
  expect_equal(e$total, 1.5, tolerance = 1e-7)
  expect_true(all(e$violations$violation > 0))
  expect_output(print(e), "must move")
  # the model itself is unchanged
  expect_length(m$cons, 3L)
  expect_equal(nrow(opt_variables(m)), 2L)
  expect_equal(opt_status(opt_solve(m)), "infeasible")
  # with a heavy row weight, the upper bound of x[1] moves instead
  m2 <- opt_model()
  x2 <- opt_variable(m2, "x", dim = 2, lower = 0, upper = c(1, 5))
  opt_constraint(m2, x2[1] >= 2, name = "demand_1")
  e2 <- opt_explain_infeasibility(m2, weights = c(demand_1 = 100))
  expect_equal(e2$violations$kind[1], "bound")
  expect_equal(e2$violations$name[1], "x[1]")
  expect_equal(e2$violations$side[1], "upper")
  expect_equal(e2$violations$bound[1], 1)
  # a feasible model needs no move
  opt_delete_constraint(m, "demand_1")
  f <- opt_explain_infeasibility(m)
  expect_true(f$feasible)
  expect_equal(nrow(f$violations), 0L)
  expect_output(print(f), "feasible")
  expect_error(opt_explain_infeasibility(m, weights = c(capacity = -1)),
               "positive")
})

test_that("infeasibility explanations work for integer and nonlinear models", {
  m <- opt_model()
  n <- opt_variable(m, "n", type = "integer", lower = 0, upper = 10)
  opt_constraint(m, 2 * n == 3, name = "odd")
  e <- opt_explain_infeasibility(m, bounds = FALSE)
  expect_false(e$feasible)
  expect_equal(e$total, 1, tolerance = 1e-7)
  skip_if_not_installed("nloptr")
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 0, upper = 2)
  opt_constraint(m, exp(x) >= 20, name = "growth")
  e <- opt_explain_infeasibility(m, bounds = FALSE)
  expect_equal(e$violations$name, "growth")
  expect_equal(e$violations$violation, 20 - exp(2), tolerance = 1e-5)
})

test_that("CVXR conversion reproduces optcompR solutions", {
  skip_if_not_installed("CVXR")
  cmp <- function(m) {
    cv <- opt_to_cvxr(m)
    v <- CVXR::psolve(cv$problem)
    r <- opt_solve(m)
    expect_equal(v, opt_objective_value(r), tolerance = 1e-5)
    expect_equal(unname(opt_cvxr_values(cv)), unname(r$x), tolerance = 1e-3)
  }
  m <- opt_model("max")
  x <- opt_variable(m, "x", dim = 2, lower = 0)
  opt_constraint(m, x[1] + 2 * x[2] <= 4)
  opt_constraint(m, 3 * x[1] + x[2] <= 6)
  opt_objective(m, x[1] + x[2] + 1)
  cmp(m)
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 2)
  opt_constraint(m, x[1] + x[2] == 1)
  opt_objective(m, x[1]^2 + x[1] * x[2] + 2 * x[2]^2 - x[1])
  cmp(m)
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 2)
  opt_constraint(m, x[1]^2 + x[2]^2 <= 1)
  opt_objective(m, -x[1] - 2 * x[2])
  cmp(m)
  m <- opt_model()
  x <- opt_variable(m, "x", dim = 3, type = "integer", lower = 0, upper = 5)
  opt_constraint(m, 2 * x[1] + 3 * x[2] + x[3] >= 7.5)
  opt_objective(m, sum(x * c(3, 4, 2)))
  cmp(m)
  m <- opt_model()
  t <- opt_variable(m, "t")
  y <- opt_variable(m, "y", dim = 2)
  opt_constraint(m, y[1] + y[2] == 2)
  opt_soc_constraint(m, y - c(3, 0), t)
  opt_objective(m, t)
  cmp(m)
  m <- opt_model()
  X <- opt_variable(m, "X", dim = c(2, 2))
  opt_psd_constraint(m, X)
  opt_constraint(m, X[1, 1] + X[2, 2] == 1)
  opt_constraint(m, X[1, 2] == X[2, 1])
  opt_objective(m, -(2 * X[1, 2] + X[1, 1]))
  cmp(m)
  m <- opt_model()
  u <- opt_variable(m, "u", dim = 3)
  opt_pow_cone(m, u[1], u[2], u[3], 0.3)
  opt_constraint(m, u[1] + u[2] <= 2)
  opt_objective(m, -u[3])
  cmp(m)
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 1)
  opt_objective(m, exp(x))
  expect_error(opt_to_cvxr(m), "general nonlinear")
  expect_error(opt_cvxr_values(list()), "opt_to_cvxr")
})

test_that("SCIP solves polynomial MINLPs globally by exact lifting", {
  skip_if_not_installed("scip")
  m <- opt_model()
  x <- opt_variable(m, "x", lower = -2, upper = 2)
  opt_objective(m, x^4 - 3 * x^2 + x)
  r <- opt_solve(m)
  expect_equal(r$solver, "scip")
  expect_equal(opt_status(r), "optimal")
  xs <- seq(-2, 2, length.out = 200001)
  expect_equal(opt_objective_value(r), min(xs^4 - 3 * xs^2 + xs),
               tolerance = 1e-6)
  expect_length(r$x, 1L)
  m <- opt_model("max")
  x <- opt_variable(m, "x", lower = 0, upper = 3)
  n <- opt_variable(m, "n", type = "integer", lower = 0, upper = 4)
  opt_constraint(m, x^3 + n * x^2 <= 10)
  opt_constraint(m, x * n * x >= 1)
  opt_objective(m, x + n)
  r <- opt_solve(m)
  expect_equal(r$solver, "scip")
  expect_equal(unname(opt_value(r, n)), 4)
  xv <- unname(opt_value(r, x))
  expect_equal(xv^3 + 4 * xv^2, 10, tolerance = 1e-6)
  expect_true(opt_check(r)$feasible)
  # powers with different exponents per element and division by constants
  m <- opt_model()
  v <- opt_variable(m, "v", dim = 2, lower = -1, upper = 2)
  opt_objective(m, sum(v^c(3, 4)) / 2 - v[1] * v[2])
  r <- opt_solve(m, solver = "scip")
  g <- expand.grid(a = seq(-1, 2, by = 0.005), b = seq(-1, 2, by = 0.005))
  expect_lte(opt_objective_value(r),
             min((g$a^3 + g$b^4) / 2 - g$a * g$b) + 1e-9)
  # non-polynomial operators are refused with the reason
  m <- opt_model()
  x <- opt_variable(m, "x", lower = -1, upper = 1)
  opt_objective(m, sin(3 * x) + x^4)
  expect_false(identical(opt_solve(m)$solver, "scip"))
  expect_error(opt_solve(m, solver = "scip"), "polynomial nonlinear models only")
  # convex nonlinear models still go to local NLP solvers automatically
  m <- opt_model()
  x <- opt_variable(m, "x", lower = 0.5, upper = 2)
  opt_objective(m, x^4 + x)
  expect_false(identical(opt_solve(m)$solver, "scip"))
  expect_equal(opt_solve(m, solver = "scip")$solver, "scip")
})
