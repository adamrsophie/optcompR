# Model building and solve times for optcompR, with ompr for comparison
# when it is installed. Run with: Rscript benchmarks/run.R
library(optcompR)
library(Matrix)

time_it <- function(expr) unname(system.time(expr)[["elapsed"]])

bench_lp <- function(n, k) {
  set.seed(n)
  A <- rsparsematrix(k, n, density = 5 / n, rand.x = function(m) runif(m, 0.5, 2)) +
    sparseMatrix(i = seq_len(k), j = sample(n, k), x = 1, dims = c(k, n))
  b <- runif(k, 1, 10)
  cost <- runif(n, 1, 3)
  m <- opt_model()
  build <- time_it({
    x <- opt_variable(m, "x", dim = n, lower = 0)
    opt_constraint(m, A %*% x >= b)
    opt_objective(m, sum(cost * x))
  })
  r <- NULL
  solve <- time_it(r <- opt_solve(m))
  data.frame(problem = "sparse LP", n = n, rows = k, build = build,
             solve = solve, status = opt_status(r))
}

bench_scalar <- function(k) {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = k, lower = 0, upper = 1)
  build <- time_it(for (i in seq_len(k - 1L)) opt_constraint(m, x[i] + x[i + 1L] <= 1))
  opt_objective(m, sum(x), sense = "max")
  r <- NULL
  solve <- time_it(r <- opt_solve(m))
  data.frame(problem = "scalar constraints", n = k, rows = k - 1L,
             build = build, solve = solve, status = opt_status(r))
}

bench_ompr <- function(k) {
  if (!requireNamespace("ompr", quietly = TRUE) ||
      !requireNamespace("ompr.roi", quietly = TRUE) ||
      !requireNamespace("ROI.plugin.glpk", quietly = TRUE)) return(NULL)
  library(ompr)
  build <- time_it({
    mo <- MILPModel() |>
      add_variable(x[i], i = 1:k, lb = 0, ub = 1) |>
      set_objective(sum_expr(x[i], i = 1:k), "max") |>
      add_constraint(x[i] + x[i + 1] <= 1, i = 1:(k - 1))
  })
  solve <- time_it(res <- solve_model(mo, ompr.roi::with_ROI("glpk")))
  data.frame(problem = "ompr scalar constraints", n = k, rows = k - 1L,
             build = build, solve = solve, status = res$status)
}

out <- rbind(bench_lp(1e4, 1e3), bench_lp(1e5, 1e4), bench_scalar(1e3),
             bench_scalar(1e4), bench_ompr(1e3), bench_ompr(1e4))
print(out, row.names = FALSE)
