# Benchmarks of optcompR: model building, solves through optcompR against
# the same solver called directly and through ROI, repeated solves after
# parameter updates, solution derivatives and memory.
#
# Run with: Rscript benchmarks/run.R [output directory]
# Results are printed and written to results.csv and results.md.
suppressPackageStartupMessages({
  library(optcompR)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[[1L]] else "benchmarks"
has <- function(...) all(vapply(c(...), requireNamespace, TRUE, quietly = TRUE))

# median elapsed seconds over `reps` runs, and memory allocated by one run
measure <- function(expr, reps = 3L) {
  # the benchmark expression, evaluated in the caller's frame
  e <- substitute(expr)
  env <- parent.frame()
  f <- function() eval(e, env)
  times <- vapply(seq_len(reps), function(i)
    unname(system.time(f())[["elapsed"]]), 0)
  mem <- if (has("bench")) as.numeric(bench::bench_memory(f())$mem_alloc) else NA
  c(seconds = stats::median(times), mb = mem / 2^20)
}

rows <- list()
add <- function(group, task, size, way, m) {
  message(sprintf("%-12s %-40s %-10s %-22s %8.3f s", group, task, size, way,
                  m[["seconds"]]))
  rows[[length(rows) + 1L]] <<- data.frame(
    group = group, task = task, size = size, way = way,
    seconds = round(m[["seconds"]], 4), mb = round(m[["mb"]], 1))
}

# ----------------------------------------------------------- sparse LP
lp_data <- function(n, k) {
  set.seed(n)
  A <- rsparsematrix(k, n, density = 5 / n,
                     rand.x = function(m) runif(m, 0.5, 2)) +
    sparseMatrix(i = seq_len(k), j = sample(n, k), x = 1, dims = c(k, n))
  list(A = A, b = runif(k, 1, 10), cost = runif(n, 1, 3))
}
lp_model <- function(d) {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = ncol(d$A), lower = 0)
  opt_constraint(m, d$A %*% x >= d$b)
  opt_objective(m, sum(d$cost * x))
  m
}
for (sz in list(c(1e4, 1e3), c(5e4, 5e3))) {
  d <- lp_data(sz[1], sz[2])
  size <- sprintf("%g x %g", sz[2], sz[1])
  add("LP", "build", size, "optcompR", measure(lp_model(d)))
  add("LP", "build and solve", size, "optcompR highs",
      measure(opt_solve(lp_model(d), solver = "highs")))
  add("LP", "build and solve", size, "highs direct", measure(
    highs::highs_solve(L = d$cost, lower = rep(0, ncol(d$A)), upper = Inf,
                       A = d$A, lhs = d$b, rhs = rep(Inf, nrow(d$A)))))
  if (has("ROI", "ROI.plugin.highs")) {
    add("LP", "build and solve", size, "ROI highs", measure({
      T <- summary(d$A)
      st <- slam::simple_triplet_matrix(T$i, T$j, T$x, nrow(d$A), ncol(d$A))
      op <- ROI::OP(d$cost, ROI::L_constraint(st, rep(">=", nrow(d$A)), d$b))
      ROI::ROI_solve(op, solver = "highs")
    }))
  }
}

# ---------------------------------------------------------- convex QP
qp_model <- function(n) {
  set.seed(n)
  F <- rsparsematrix(n, n, density = 3 / n)
  P <- crossprod(F) + Diagonal(n, 0.1)
  q <- rnorm(n)
  m <- opt_model()
  x <- opt_variable(m, "x", dim = n, lower = -1, upper = 1)
  opt_constraint(m, sum(x) == 0)
  opt_objective(m, 0.5 * sum(x * (P %*% x)) + sum(q * x))
  list(m = m, P = P, q = q)
}
for (n in c(1e3, 5e3)) {
  qp <- qp_model(n)
  size <- format(n, scientific = FALSE)
  add("QP", "build", size, "optcompR", measure(qp_model(n)))
  for (s in intersect(c("highs", "osqp", "clarabel"),
                      with(opt_solvers(), solver[available]))) {
    add("QP", "solve", size, paste("optcompR", s),
        measure(opt_solve(qp$m, solver = s)))
  }
  if (has("osqp")) {
    A <- rbind(Matrix(1, 1, n, sparse = TRUE), Diagonal(n))
    add("QP", "solve", size, "osqp direct", measure(osqp::solve_osqp(
      qp$P, qp$q, A, c(0, rep(-1, n)), c(0, rep(1, n)),
      osqp::osqpSettings(verbose = FALSE))))
  }
}

# ------------------------------------------------------- knapsack MILP
for (n in c(50, 200)) {
  set.seed(n)
  w <- runif(n, 1, 10)
  v <- w + runif(n, 0, 5)
  cap <- sum(w) / 3
  size <- format(n)
  build <- function() {
    m <- opt_model("max")
    x <- opt_variable(m, "x", dim = n, type = "binary")
    opt_constraint(m, sum(w * x) <= cap)
    opt_objective(m, sum(v * x))
    m
  }
  for (s in intersect(c("highs", "scip"), with(opt_solvers(), solver[available]))) {
    add("MILP", "build and solve", size, paste("optcompR", s),
        measure(opt_solve(build(), solver = s)))
  }
  add("MILP", "build and solve", size, "highs direct", measure(
    highs::highs_solve(L = -v, lower = rep(0, n), upper = rep(1, n),
                       A = matrix(w, 1), lhs = -Inf, rhs = cap,
                       types = rep(2L, n))))
  if (has("ROI", "ROI.plugin.glpk")) {
    add("MILP", "build and solve", size, "ROI glpk", measure({
      op <- ROI::OP(v, ROI::L_constraint(matrix(w, 1), "<=", cap),
                    types = rep("B", n), maximum = TRUE)
      ROI::ROI_solve(op, solver = "glpk")
    }))
  }
}

# -------------------------------------------------------------- NLP
nlp_model <- function(n) {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = n, lower = -5, upper = 5, start = 0.5)
  opt_constraint(m, x[-1] - x[-n]^2 / 4 == 0)
  opt_objective(m, sum((x - 1)^2) + sum(exp(x / 10)))
  m
}
for (n in c(50, 400)) {
  size <- format(n)
  for (s in intersect(c("ipopt", "nloptr"), with(opt_solvers(), solver[available]))) {
    if (s == "nloptr" && n > 100) next # dense quasi-Newton, too slow
    add("NLP", "build and solve", size, paste("optcompR", s),
        measure(opt_solve(nlp_model(n), solver = s), reps = 1L))
  }
}

# ----------------------------------------------------- repeated solves
rep_model <- function(n) {
  set.seed(1)
  A <- matrix(runif(40 * n), 40)
  m <- opt_model("max")
  x <- opt_variable(m, "x", dim = n, lower = 0, upper = 1)
  b <- opt_parameter(m, "b", dim = 40, value = rep(10, 40))
  opt_constraint(m, A %*% x <= b)
  opt_objective(m, sum(runif(n) * x))
  list(m = m, b = b, A = A)
}
rm_ <- rep_model(500)
add("repeated", "20 solves after parameter updates", "40 x 500",
    "optcompR parameters", measure({
      for (k in 1:20) {
        opt_set_parameter(rm_$m, rm_$b, rep(10 + k / 10, 40))
        opt_solve(rm_$m, solver = "highs")
      }
    }))
add("repeated", "20 solves after parameter updates", "40 x 500",
    "rebuild each time", measure({
      for (k in 1:20) {
        r <- rep_model(500)
        opt_set_parameter(r$m, r$b, rep(10 + k / 10, 40))
        opt_solve(r$m, solver = "highs")
      }
    }))
scen <- lapply(1:20, function(k) rep(10 + k / 10, 40))
add("repeated", "20 scenarios", "40 x 500", "opt_scenarios",
    measure(opt_scenarios(rm_$m, rm_$b, scen, solver = "highs")))

# --------------------------------------------------------- derivatives
der_model <- function(n) {
  set.seed(2)
  m <- opt_model()
  x <- opt_variable(m, "x", dim = n)
  p <- opt_parameter(m, "p", dim = n, value = rnorm(n))
  opt_constraint(m, sum(x) == 1)
  opt_objective(m, sum((x - p)^2) + 0.1 * sum(x^4))
  list(m = m, p = p)
}
for (n in c(20, 200)) {
  dm <- der_model(n)
  r <- opt_solve(dm$m)
  add("derivatives", "solution Jacobian", format(n), "implicit KKT",
      measure(opt_solution_jacobian(r, dm$p)))
}

# ------------------------------------------------- scalar construction
scalar_build <- function(k) {
  m <- opt_model()
  x <- opt_variable(m, "x", dim = k, lower = 0, upper = 1)
  for (i in seq_len(k - 1L)) opt_constraint(m, x[i] + x[i + 1L] <= 1)
  m
}
for (k in c(1e3, 5e3)) {
  add("construction", "scalar constraints in a loop", format(k), "optcompR",
      measure(scalar_build(k), reps = 1L))
  add("construction", "one vector constraint", format(k), "optcompR",
      measure({
        m <- opt_model()
        x <- opt_variable(m, "x", dim = k, lower = 0, upper = 1)
        opt_constraint(m, x[-k] + x[-1] <= 1)
      }))
  # ompr loads its own OpenMP runtime on some systems, which aborts R when
  # another one is loaded, so the comparison is opt-in
  if (k <= 1e3 && identical(Sys.getenv("OPTCOMPR_BENCH_OMPR"), "true") &&
      has("ompr", "ompr.roi")) {
    add("construction", "scalar constraints in a loop", format(k), "ompr",
        measure({
          ompr::MILPModel() |>
            ompr::add_variable(x[i], i = 1:k, lb = 0, ub = 1) |>
            ompr::add_constraint(x[i] + x[i + 1] <= 1, i = 1:(k - 1))
        }, reps = 1L))
  }
}

# scalar variables, quadratic rows and indicators, one at a time or as
# blocks: paths whose cost must stay linear in the model size
for (k in c(2e3, 8e3)) {
  add("construction", "scalar variables in a loop", format(k), "optcompR",
      measure({
        m <- opt_model()
        for (i in seq_len(k)) opt_variable(m, paste0("v", i))
      }, reps = 1L))
  add("construction", "quadratic rows x_i^2 <= 1", format(k), "optcompR",
      measure({
        m <- opt_model()
        x <- opt_variable(m, "x", dim = k)
        opt_constraint(m, x^2 <= 1)
        opt_standard_form(m)
      }))
  add("construction", "indicators z_i => x_i <= 1", format(k), "optcompR",
      measure({
        m <- opt_model("max")
        x <- opt_variable(m, "x", dim = k, lower = 0, upper = 10)
        z <- opt_variable(m, "z", dim = k, type = "binary")
        opt_indicator(m, z, x <= 1)
        opt_objective(m, sum(x) - sum(z))
        opt_standard_form(m)
      }))
}

# repeated convex QP solves: the OSQP workspace is kept and updated
if (has("osqp")) {
  n <- 2000
  set.seed(3)
  m <- opt_model()
  x <- opt_variable(m, "x", dim = n, lower = -1, upper = 1)
  a <- opt_parameter(m, "a", dim = n, value = rnorm(n))
  opt_objective(m, sum((x - a)^2))
  vals <- lapply(1:20, function(i) rnorm(n))
  add("repeated", "20 QP solves after parameter updates", format(n),
      "optcompR osqp kept", measure({
        for (v in vals) {
          opt_set_parameter(m, a, v)
          opt_solve(m, solver = "osqp")
        }
      }))
  add("repeated", "20 QP solves after parameter updates", format(n),
      "optcompR osqp new", measure({
        for (v in vals) {
          opt_set_parameter(m, a, v)
          opt_solve(m, solver = "osqp", warm_start = FALSE)
        }
      }))
}

out <- do.call(rbind, rows)
print(out, row.names = FALSE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
utils::write.csv(out, file.path(out_dir, "results.csv"), row.names = FALSE)
md <- c(sprintf("R %s, %s, optcompR %s.", getRversion(), R.version$platform,
                packageVersion("optcompR")), "",
        "| group | task | size | way | seconds | MB |",
        "|---|---|---|---|---:|---:|",
        sprintf("| %s | %s | %s | %s | %s | %s |", out$group, out$task,
                out$size, out$way, format(out$seconds), format(out$mb)))
writeLines(md, file.path(out_dir, "results.md"))
