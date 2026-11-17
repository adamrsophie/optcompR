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
  opts <- control$options
  ctrl_args <- utils::modifyList(
    list(threads = 1L, log_to_console = isTRUE(control$verbose)), opts)
  # HiGHS keeps one global thread scheduler per process and cannot resize
  # it: a different thread count aborts R on Windows even after
  # hi_reset_global_scheduler(), so the first count is kept for the session.
  threads <- as.integer(ctrl_args$threads)
  first <- .optcompR_state$highs_threads
  if (is.null(first)) {
    .optcompR_state$highs_threads <- threads
  } else if (!identical(first, threads)) {
    stop(sprintf(paste(
      "HiGHS already runs with %d thread(s) in this R session and cannot",
      "change its thread count; restart R to use threads = %d."),
      first, threads), call. = FALSE)
  }
  if (is.finite(control$time_limit)) ctrl_args$time_limit <- control$time_limit
  if (!is.null(control$mip_gap)) ctrl_args$mip_rel_gap <- control$mip_gap
  # the native solver is kept per model and reused while the structure,
  # the sparsity pattern, Q and the controls stay the same: costs, bounds
  # and matrix values are then updated in place, and HiGHS restarts from
  # its last basis
  cache <- control$cache
  data <- list(c = s * sf$c, lower = sf$lower, upper = sf$upper,
               lhs = lhs %||% numeric(), rhs = rhs %||% numeric(),
               ax = if (is.null(A)) numeric() else A@x, offset = s * sf$offset)
  sig <- list(key = sf$structure_key, n = n, m = nrow(sf$A), types = types,
              ai = if (is.null(A)) NULL else A@i,
              ap = if (is.null(A)) NULL else A@p, Q = Q, ctrl = ctrl_args)
  sess <- if (isTRUE(control$warm_start) && is.environment(cache)) cache$highs
  reused <- !is.null(sess) && identical(sess$sig, sig) &&
    highs_update(sess, data)
  if (reused) {
    solver <- sess$solver
  } else {
    hm <- highs::highs_model(Q = Q, L = data$c, lower = sf$lower,
                             upper = sf$upper, A = A, lhs = lhs, rhs = rhs,
                             types = types, maximum = FALSE,
                             offset = data$offset)
    ctrl <- do.call(highs::highs_control, ctrl_args)
    solver <- highs::highs_solver(hm, control = ctrl)
  }
  if (is.environment(cache)) {
    cache$highs <- list(solver = solver, sig = sig, data = data)
  }
  start_used <- FALSE
  if (isTRUE(control$warm_start) && any(!is.na(sf$start))) {
    # partial starts are allowed: NA entries are left to HiGHS
    highs::hi_solver_set_start(solver$solver, sf$start)
    start_used <- TRUE
  }
  # an LP re-solved after a parameter update starts from the previous
  # optimal basis when the structure (and so the basis size) is unchanged
  cache <- control$cache
  bkey <- paste(sf$structure_key, n, nrow(sf$A))
  basis_used <- FALSE
  if (!reused && !cls$integer && is.null(sf$Q) && isTRUE(control$warm_start) &&
      !start_used && is.environment(cache) &&
      identical(cache$highs_basis_key, bkey)) {
    b <- cache$highs_basis
    ok <- tryCatch({
      highs::hi_solver_set_basis(solver$solver, b$col_status, b$row_status)
      TRUE
    }, error = function(e) FALSE)
    basis_used <- isTRUE(ok)
  }
  solver$solve()
  code <- solver$status()
  if (!cls$integer && is.null(sf$Q) && is.environment(cache) &&
      highs_termination(code) == "optimal") {
    b <- tryCatch(highs::hi_solver_get_basis(solver$solver),
                  error = function(e) NULL)
    if (isTRUE(b$valid)) {
      cache$highs_basis <- b
      cache$highs_basis_key <- bkey
    }
  }
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
  # ranging is computed on request by opt_sensitivity(), from the basis
  basis <- NULL
  if (cls$class == "LP" && term == "optimal") {
    basis <- tryCatch(highs::hi_solver_get_basis(solver$solver),
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
       local = FALSE, warm_start = start_used || basis_used || reused,
       extra = list(basis = basis, sign = s, reused = reused),
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

# Apply the changed costs, bounds, matrix values and offset to a kept HiGHS
# solver. FALSE when the solver can no longer be used (for example after the
# model was serialized), so a new one is built.
highs_update <- function(sess, data) {
  old <- sess$data
  sv <- sess$solver
  tryCatch({
    highs::hi_solver_status(sv$solver)
    d <- which(old$c != data$c)
    if (length(d)) sv$L(d, data$c[d])
    d <- which(old$lower != data$lower | old$upper != data$upper)
    if (length(d)) sv$vbounds(d, data$lower[d], data$upper[d])
    d <- which(old$lhs != data$lhs | old$rhs != data$rhs)
    if (length(d)) sv$cbounds(d, data$lhs[d], data$rhs[d])
    d <- which(old$ax != data$ax)
    if (length(d)) {
      # the column of entry k of a column-compressed matrix with pointers p
      p <- sess$sig$ap
      jj <- findInterval(d - 1L, p, left.open = FALSE)
      sv$A(sess$sig$ai[d] + 1L, jj, data$ax[d])
    }
    if (!identical(old$offset, data$offset)) {
      highs::hi_solver_set_offset(sv$solver, data$offset)
    }
    TRUE
  }, error = function(e) FALSE)
}

# Sensitivity ranging of an optimal HiGHS LP result, from a solver rebuilt
# from the result's standard form and started at its optimal basis.
highs_ranging <- function(result) {
  sf <- result$standard_form
  b <- result$extra$basis
  if (is.null(b)) return(NULL)
  s <- result$extra$sign
  A <- if (nrow(sf$A)) methods::as(sf$A, "CsparseMatrix") else NULL
  hm <- highs::highs_model(L = s * sf$c, lower = sf$lower, upper = sf$upper,
                           A = A, lhs = if (!is.null(A)) sf$row_lower,
                           rhs = if (!is.null(A)) sf$row_upper,
                           maximum = FALSE, offset = s * sf$offset)
  args <- list(threads = .optcompR_state$highs_threads %||% 1L,
               log_to_console = FALSE)
  solver <- highs::highs_solver(hm, control = do.call(highs::highs_control, args))
  tryCatch(highs::hi_solver_set_basis(solver$solver, b$col_status, b$row_status),
           error = function(e) NULL)
  solver$solve()
  tryCatch(highs::hi_solver_get_ranging(solver$solver),
           error = function(e) NULL)
}
