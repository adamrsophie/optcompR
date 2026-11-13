# HiGHS backend through the 'highs' package.
#
# Models are always passed in minimization form: the objective is multiplied
# by s = -1 for maximization and the results are multiplied back. HiGHS row
# and column duals of a minimization problem equal the derivative of the
# optimal objective with respect to the active bound (checked by finite
# differences in the tests), so multiplying by s gives the shadow price
# d f* / d b of the user's problem in both senses.

highs_termination <- function(code) {
  switch(as.character(code),
    "7" = "optimal", "8" = "infeasible", "9" = "infeasible_or_unbounded",
    "10" = "unbounded", "11" = "objective_limit", "12" = "objective_limit",
    "13" = "time_limit", "14" = "iteration_limit", "16" = "solution_limit",
    "17" = "interrupted", "18" = "memory_limit",
    "1" = , "2" = "invalid_model", "3" = , "4" = , "5" = "numerical_error",
    "unknown")
}

highs_point_status <- function(s) {
  switch(s %||% "", Feasible = "feasible_point", Infeasible = "infeasible_point",
         None = "no_solution", "unknown")
}

highs_solve_sf <- function(sf, cls, control) {
  s <- if (sf$sense == "max") -1 else 1
  n <- length(sf$c)
  types <- ifelse(sf$type %in% c("integer", "binary"), 2L, 1L)
  A <- NULL
  lhs <- NULL
  rhs <- NULL
  if (nrow(sf$A) > 0L) {
    A <- methods::as(sf$A, "CsparseMatrix")
    lhs <- sf$row_lower
    rhs <- sf$row_upper
  }
  Q <- if (is.null(sf$Q)) NULL else
    methods::as(methods::as(s * sf$Q, "generalMatrix"), "CsparseMatrix")
  hm <- highs::highs_model(Q = Q, L = s * sf$c, lower = sf$lower,
                           upper = sf$upper, A = A, lhs = lhs, rhs = rhs,
                           types = types, maximum = FALSE,
                           offset = s * sf$offset)
  opts <- control$options
  ctrl_args <- c(list(threads = 1L, log_to_console = isTRUE(control$verbose)),
                 opts)
  if (is.finite(control$time_limit)) ctrl_args$time_limit <- control$time_limit
  if (!is.null(control$mip_gap)) ctrl_args$mip_rel_gap <- control$mip_gap
  ctrl <- do.call(highs::highs_control, ctrl_args)
  solver <- highs::highs_solver(hm, control = ctrl)
  start_used <- FALSE
  if (isTRUE(control$warm_start) && any(!is.na(sf$start))) {
    st <- sf$start
    st[is.na(st)] <- pmin(pmax(0, sf$lower[is.na(st)]), sf$upper[is.na(st)])
    ok <- tryCatch({
      highs::hi_solver_set_start(solver$solver, list(col_value = st))
      TRUE
    }, error = function(e) FALSE)
    start_used <- ok
  }
  solver$solve()
  code <- solver$status()
  sol <- solver$solution()
  info <- solver$info()
  integer <- cls$integer
  term <- highs_termination(code)
  primal <- highs_point_status(info$primal_solution_status)
  if (!isTRUE(sol$value_valid)) primal <- "no_solution"
  dual <- if (integer) "not_available" else
    highs_point_status(info$dual_solution_status)
  x <- if (primal == "no_solution") rep(NA_real_, n) else sol$col_value
  rd <- if (!integer && isTRUE(sol$dual_valid)) s * sol$row_dual else NULL
  cd <- if (!integer && isTRUE(sol$dual_valid)) s * sol$col_dual else NULL
  ranging <- NULL
  if (cls$class == "LP" && term == "optimal") {
    ranging <- tryCatch(highs::hi_solver_get_ranging(solver$solver),
                        error = function(e) NULL)
  }
  obj <- if (primal == "no_solution") NA_real_ else
    s * info$objective_function_value
  gap <- if (integer) info$mip_gap else NA_real_
  bound <- if (integer && is.finite(info$mip_dual_bound)) {
    s * info$mip_dual_bound
  } else if (!integer && term == "optimal") obj else NA_real_
  list(termination = term, primal_status = primal, dual_status = dual,
       x = x, objective = obj, row_dual = rd, col_dual = cd,
       iterations = c(simplex = info$simplex_iteration_count,
                      ipm = info$ipm_iteration_count,
                      qp = info$qp_iteration_count,
                      nodes = info$mip_node_count),
       gap = gap, bound = bound, message = solver$status_message(),
       local = FALSE, warm_start = start_used,
       extra = list(ranging = ranging, sign = s),
       raw = list(status = code, solution = sol, info = info))
}

highs_version <- function() {
  as.character(utils::packageVersion("highs"))
}

register_highs <- function() {
  register_solver(
    "highs", "highs",
    "HiGHS: dual simplex, interior point and branch and bound (LP, convex QP, MILP)",
    classes = c("LP", "MILP", "QP"), convex_only = "QP",
    global = c("LP", "MILP", "QP"),
    features = list(duals = TRUE, ranging = TRUE, warm_start = TRUE,
                    time_limit = TRUE, mip_gap = TRUE, iteration_limit = TRUE),
    solve = highs_solve_sf, version = highs_version, priority = 10L)
}
