# Optional interface to the Z3 SMT solver through z3r.
#
# Only linear models are translated: linear integer and real arithmetic is
# decidable, so "sat" and "unsat" answers are exact. Nonlinear arithmetic in
# Z3 is incomplete and has different semantics from numerical optimization,
# so quadratic and nonlinear models are refused.

need_z3r <- function() {
  if (!requireNamespace("z3r", quietly = TRUE)) {
    stop("This function needs the z3r package (remotes::install_github('adamrsophie/z3r')).",
         call. = FALSE)
  }
}

# Exact rational text for a coefficient.
z3_number_text <- function(v, coefficients) {
  if (coefficients == "decimal") {
    return(format(v, digits = 15L, scientific = FALSE, trim = TRUE))
  }
  # binary: the exact value of the double, num / 2^k
  k <- 0L
  while (v != round(v) && k < 1100L) {
    v <- v * 2
    k <- k + 1L
  }
  if (k == 0L) sprintf("%.0f", v) else sprintf("%.0f/%.0f", v, 2^k)
}

z3_translate <- function(sf, coefficients) {
  if (any(sf$row_kind != "linear")) {
    stop(paste("SMT translation supports linear constraints only; Z3",
               "nonlinear arithmetic is incomplete."), call. = FALSE)
  }
  ctx <- z3r::z3_context()
  n <- length(sf$c)
  vars <- lapply(seq_len(n), function(j) {
    nm <- sprintf("v%d", j)
    if (sf$type[j] %in% c("integer", "binary")) z3r::z3_int(ctx, nm) else
      z3r::z3_real(ctx, nm)
  })
  num <- function(v) z3r::z3_real_val(ctx, z3_number_text(v, coefficients))
  facts <- list()
  for (j in seq_len(n)) {
    if (is.finite(sf$lower[j])) facts <- c(facts, list(vars[[j]] >= num(sf$lower[j])))
    if (is.finite(sf$upper[j])) facts <- c(facts, list(vars[[j]] <= num(sf$upper[j])))
  }
  lin <- function(cols, coefs, const = 0) {
    terms <- Map(function(j, a) num(a) * vars[[j]], cols, coefs)
    if (const != 0 || length(terms) == 0L) terms <- c(terms, list(num(const)))
    if (length(terms) == 1L) terms[[1L]] else do.call(z3r::z3_sum, terms)
  }
  rt <- row_terms(sf$A)
  for (r in seq_len(nrow(sf$A))) {
    idx <- if (rt$p[r + 1L] > rt$p[r]) (rt$p[r] + 1L):rt$p[r + 1L] else integer()
    e <- lin(rt$j[idx], rt$x[idx])
    lo <- sf$row_lower[r]
    hi <- sf$row_upper[r]
    if (is.finite(lo) && is.finite(hi) && lo == hi) {
      facts <- c(facts, list(e == num(lo)))
    } else {
      if (is.finite(lo)) facts <- c(facts, list(e >= num(lo)))
      if (is.finite(hi)) facts <- c(facts, list(e <= num(hi)))
    }
  }
  nz <- which(sf$c != 0)
  list(ctx = ctx, vars = vars, facts = facts,
       objective = lin(nz, sf$c[nz], sf$offset))
}

z3_run <- function(tr, extra = list()) {
  s <- z3r::z3_solver(tr$ctx)
  for (f in c(tr$facts, extra)) z3r::z3_add(s, f)
  st <- z3r::z3_check(s)
  vals <- NULL
  if (identical(st, "sat")) {
    mo <- z3r::z3_model(s)
    vals <- vapply(tr$vars, function(v) {
      as.numeric(z3r::z3_value(z3r::z3_eval(mo, v, completion = TRUE),
                               exact = FALSE))
    }, 0)
  }
  list(status = st, values = vals)
}

#' Exact feasibility check with Z3
#'
#' Translates the constraints and bounds of a linear model (LP or MILP,
#' after parameter substitution) to linear integer and real arithmetic and
#' asks Z3 whether a feasible point exists. For these theories Z3 is a
#' decision procedure, so `"sat"` and `"unsat"` are exact answers about the
#' model as written, without floating-point tolerances. A satisfying
#' assignment is a feasible point, not an optimal one.
#'
#' Coefficients and bounds are doubles. With `coefficients = "decimal"` each
#' is read as its 15-significant-digit decimal (so `0.1` means 1/10); with
#' `"binary"` as the exact binary value of the double.
#'
#' @param model An [opt_model()] with linear constraints.
#' @param objective_bound Optional number: also require the objective to be
#'   at most (minimization) or at least (maximization) this value.
#' @param coefficients `"decimal"` or `"binary"`, see Details.
#' @return A list with `status` (`"sat"`, `"unsat"` or `"unknown"`) and, when
#'   satisfiable, `values`, a named feasible point.
#' @seealso [opt_z3_certify()]
#' @examplesIf requireNamespace("z3r", quietly = TRUE)
#' m <- opt_model()
#' x <- opt_variable(m, "x", type = "integer", lower = 0, upper = 10)
#' opt_constraint(m, 2 * x == 7)
#' opt_z3_check(m)$status
#' @export
opt_z3_check <- function(model, objective_bound = NULL,
                         coefficients = c("decimal", "binary")) {
  check_model(model)
  need_z3r()
  coefficients <- match.arg(coefficients)
  sf <- opt_standard_form(model)
  tr <- z3_translate(sf, coefficients)
  extra <- list()
  if (!is.null(objective_bound)) {
    if (sf$objective_kind != "linear") {
      stop("An objective bound needs a linear objective.", call. = FALSE)
    }
    b <- z3r::z3_real_val(tr$ctx, z3_number_text(objective_bound, coefficients))
    extra <- list(if (sf$sense == "min") tr$objective <= b else tr$objective >= b)
  }
  out <- z3_run(tr, extra)
  if (!is.null(out$values)) names(out$values) <- sf$col_names
  out
}

#' Certify a MILP or LP optimum with Z3
#'
#' Checks with exact arithmetic that no feasible point has an objective
#' better than the reported optimum by more than `tol * max(1, |f*|)`. An
#' `"unsat"` answer is a proof (for the model as translated, see
#' [opt_z3_check()]) that the solver's objective value is optimal up to that
#' margin. It does not prove that the returned point itself is feasible in
#' exact arithmetic; floating-point solutions satisfy constraints within
#' solver tolerances. Exact solving can be slow for large integer models.
#'
#' @param result An `opt_result` of an LP or MILP with a solution.
#' @param tol Relative margin.
#' @param coefficients See [opt_z3_check()].
#' @return A list with `certified` (`TRUE`, `FALSE` or `NA` when Z3 returns
#'   unknown), `status`, the `threshold` that was tested and, when a better
#'   point exists, `better_point`.
#' @examplesIf requireNamespace("z3r", quietly = TRUE)
#' m <- opt_model("max")
#' x <- opt_variable(m, "x", type = "binary")
#' y <- opt_variable(m, "y", type = "integer", lower = 0, upper = 3)
#' opt_constraint(m, x + y <= 4)
#' opt_objective(m, 3 * x + 2 * y)
#' opt_z3_certify(opt_solve(m))$certified
#' @export
opt_z3_certify <- function(result, tol = 1e-9,
                           coefficients = c("decimal", "binary")) {
  check_result(result)
  need_z3r()
  coefficients <- match.arg(coefficients)
  if (!result$classification$class %in% c("LP", "MILP")) {
    stop("Certification is available for LP and MILP results.", call. = FALSE)
  }
  if (!opt_has_solution(result)) stop("The result has no solution to certify.",
                                      call. = FALSE)
  sf <- result$standard_form
  v <- result$objective_value
  margin <- tol * max(1, abs(v))
  thr <- if (sf$sense == "min") v - margin else v + margin
  tr <- z3_translate(sf, coefficients)
  b <- z3r::z3_real_val(tr$ctx, z3_number_text(thr, coefficients))
  better <- if (sf$sense == "min") tr$objective < b else tr$objective > b
  out <- z3_run(tr, list(better))
  list(certified = switch(out$status, unsat = TRUE, sat = FALSE, NA),
       status = out$status, threshold = thr,
       better_point = if (!is.null(out$values))
         stats::setNames(out$values, sf$col_names))
}
