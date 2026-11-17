#' Convert a model to a CVXR problem
#'
#' Builds a CVXR `Problem` from the lowered standard form of a model
#' (parameters substituted): one CVXR variable vector for the model's
#' variables, linear and quadratic rows, variable bounds, integrality, and
#' second order, exponential, power and semidefinite cones. CVXR then checks
#' the problem against its disciplined convex programming rules, so a
#' nonconvex quadratic is refused there. Nonlinear expressions other than
#' quadratics and cones, indicators and SOS constraints have no CVXR
#' counterpart and give an error.
#'
#' @param model An [opt_model()].
#' @return An `opt_cvxr` list with `problem` (the CVXR `Problem`), `x` (the
#'   CVXR variable, one entry per model variable in the order of
#'   `col_names`), `col_names` and `offset` (the objective constant, which
#'   the CVXR objective includes).
#' @seealso [opt_cvxr_values()] to read the solution back.
#' @examplesIf requireNamespace("CVXR", quietly = TRUE)
#' m <- opt_model()
#' x <- opt_variable(m, "x", dim = 2, lower = 0)
#' opt_constraint(m, x[1] + x[2] >= 1)
#' opt_objective(m, x[1]^2 + x[2]^2)
#' cv <- opt_to_cvxr(m)
#' CVXR::psolve(cv$problem)
#' opt_cvxr_values(cv)
#' @export
opt_to_cvxr <- function(model) {
  check_model(model)
  if (!requireNamespace("CVXR", quietly = TRUE)) {
    stop("This function needs the CVXR package.", call. = FALSE)
  }
  sf <- opt_standard_form(model)
  if (identical(sf$objective_kind, "nonlinear") ||
      any(sf$row_kind == "nonlinear")) {
    stop(paste("CVXR conversion supports linear, quadratic and cone",
               "constraints only; the model has general nonlinear",
               "expressions."), call. = FALSE)
  }
  if (length(sf$indicators) || length(sf$sos)) {
    stop("Indicator and SOS constraints have no CVXR counterpart.",
         call. = FALSE)
  }
  n <- length(sf$c)
  x <- CVXR::Variable(n, name = "x")
  cons <- list()
  for (tp in c("integer", "binary")) {
    idx <- which(sf$type == tp)
    if (length(idx)) {
      z <- if (tp == "integer") CVXR::Variable(length(idx), integer = TRUE) else
        CVXR::Variable(length(idx), boolean = TRUE)
      cons <- c(cons, list(x[idx] == z))
    }
  }
  lb <- which(is.finite(sf$lower))
  ub <- which(is.finite(sf$upper))
  if (length(lb)) cons <- c(cons, list(x[lb] >= sf$lower[lb]))
  if (length(ub)) cons <- c(cons, list(x[ub] <= sf$upper[ub]))
  quad <- function(Q) 0.5 * CVXR::quad_form(x, as.matrix(Q))
  qg <- row_quad_groups(sf)
  qrows <- as.integer(names(qg))
  lin <- setdiff(seq_len(nrow(sf$A)), qrows)
  A <- as.matrix(sf$A)
  for (side in c("upper", "lower")) {
    b <- if (side == "upper") sf$row_upper else sf$row_lower
    r <- lin[is.finite(b[lin])]
    if (!length(r)) next
    eq <- r[sf$row_lower[r] == sf$row_upper[r]]
    r <- setdiff(r, eq)
    if (side == "upper" && length(eq)) {
      cons <- c(cons, list(A[eq, , drop = FALSE] %*% x == b[eq]))
    }
    if (length(r)) {
      lhs <- A[r, , drop = FALSE] %*% x
      cons <- c(cons, list(if (side == "upper") lhs <= b[r] else lhs >= b[r]))
    }
  }
  for (k in seq_along(qrows)) {
    ri <- qrows[k]
    g <- A[ri, , drop = FALSE] %*% x + quad(row_q(sf, ri, qg[[k]]))
    if (is.finite(sf$row_upper[ri])) cons <- c(cons, list(g <= sf$row_upper[ri]))
    if (is.finite(sf$row_lower[ri])) cons <- c(cons, list(g >= sf$row_lower[ri]))
  }
  for (cn in sf$cones) {
    v <- as.matrix(cn$M) %*% x + cn$const
    cons <- c(cons, list(switch(cn$type,
      soc = CVXR::SOC(v[1L], v[-1L]),
      exp = CVXR::ExpCone(v[1L], v[2L], v[3L]),
      pow = CVXR::PowCone3D(v[1L], v[2L], v[3L], cn$alpha),
      psd = {
        X <- CVXR::Variable(c(cn$k, cn$k), symmetric = TRUE)
        up <- which(upper.tri(diag(cn$k), diag = TRUE))
        list(CVXR::vec(X)[up] == v, CVXR::PSD(X))
      },
      stop(sprintf("Cone type '%s' has no CVXR counterpart.", cn$type),
           call. = FALSE))))
  }
  cons <- unlist(cons, recursive = FALSE) %||% list()
  obj <- sum(sf$c * x) + sf$offset
  if (!is.null(sf$Q) && length(sf$Q@x) && any(sf$Q@x != 0)) {
    obj <- obj + quad(sf$Q)
  }
  goal <- if (sf$sense == "max") CVXR::Maximize(obj) else CVXR::Minimize(obj)
  structure(list(problem = CVXR::Problem(goal, cons), x = x,
                 col_names = sf$col_names, offset = sf$offset),
            class = "opt_cvxr")
}

#' Read a CVXR solution back by variable name
#'
#' @param cv An `opt_cvxr` object from [opt_to_cvxr()] whose problem has
#'   been solved with `CVXR::psolve()`.
#' @return A named numeric vector, one value per model variable.
#' @export
opt_cvxr_values <- function(cv) {
  if (!inherits(cv, "opt_cvxr")) stop("`cv` must come from opt_to_cvxr().",
                                      call. = FALSE)
  v <- as.numeric(CVXR::value(cv$x))
  if (!length(v)) stop("The CVXR problem has not been solved.", call. = FALSE)
  stats::setNames(v, cv$col_names)
}

#' @export
print.opt_cvxr <- function(x, ...) {
  cat(sprintf("<opt_cvxr: %d variables>\n", length(x$col_names)))
  invisible(x)
}
