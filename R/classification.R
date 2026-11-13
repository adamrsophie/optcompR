# Problem classification with recorded evidence.

# Curvature of a symmetric matrix: is it positive semidefinite?
# Returns list(psd = TRUE/FALSE/NA, method, detail).
psd_check <- function(Q, tol = 1e-9) {
  if (is.null(Q) || length(Q@x) == 0L || all(Q@x == 0)) {
    return(list(psd = TRUE, method = "structural", min_eigen = 0,
                detail = "zero matrix"))
  }
  n <- nrow(Q)
  Qs <- Matrix::forceSymmetric(Q)
  Tq <- methods::as(Qs, "TsparseMatrix")
  if (all(Tq@i == Tq@j)) {
    d <- Matrix::diag(Q)
    return(list(psd = all(d >= 0), method = "structural", min_eigen = min(d),
                detail = "diagonal matrix; signs of the diagonal decide exactly"))
  }
  if (n <= 2000L) {
    ev <- eigen(as.matrix(Qs), symmetric = TRUE, only.values = TRUE)$values
    scale <- max(1, max(abs(ev)))
    lo <- min(ev)
    status <- if (lo >= -tol * scale) TRUE else FALSE
    detail <- sprintf("smallest eigenvalue %s (tolerance %s)", fmt_num(lo),
                      fmt_num(tol * scale))
    if (status && lo < tol * scale) {
      detail <- paste(detail, "; at the tolerance boundary, so curvature is",
                      "zero up to rounding")
    }
    return(list(psd = status, method = "numerical eigenvalues",
                min_eigen = lo, detail = detail))
  }
  shift <- tol * max(1, max(abs(Q@x)))
  ok <- tryCatch({
    Matrix::Cholesky(Qs + Matrix::Diagonal(n, shift), LDL = FALSE)
    TRUE
  }, error = function(e) FALSE)
  list(psd = if (ok) TRUE else NA, method = "numerical Cholesky",
       min_eigen = NA_real_,
       detail = if (ok) "sparse Cholesky of Q + shift I succeeded" else
         "sparse Cholesky failed; curvature undetermined")
}

#' Classify an optimization model
#'
#' Determines the problem class from the lowered standard form (after
#' parameter substitution) and records how convexity was established. Four
#' levels of evidence are distinguished: *syntactic* recognition of the
#' class, *structural* proofs (for example a diagonal quadratic form with
#' nonnegative entries), *numerical* evidence (eigenvalues compared with a
#' tolerance, which is not a formal proof near the boundary), and properties
#' that cannot be established automatically (general nonlinear expressions
#' are reported as `"unknown"`).
#'
#' Classes are `"LP"`, `"MILP"`, `"QP"`, `"MIQP"`, `"QCQP"`, `"MIQCQP"`,
#' `"NLP"` and `"MINLP"`. A model without an objective is a feasibility
#' problem and is classified by its constraints.
#'
#' @param model An [opt_model()] or an `opt_standard_form`.
#' @param tol Relative eigenvalue tolerance for curvature checks.
#' @return An `opt_classification` list with `class`, `integer`,
#'   `objective`, `constraints`, `convexity` (`"convex"`, `"nonconvex"` or
#'   `"unknown"`; for mixed-integer models this refers to the continuous
#'   relaxation) and `evidence`, a character vector.
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x")
#' y <- opt_variable(m, "y")
#' opt_objective(m, x^2 + x * y + y^2)
#' opt_classify(m)
#' opt_objective(m, x^2 - 3 * x * y + y^2)
#' opt_classify(m)$convexity
#' @export
opt_classify <- function(model, tol = 1e-9) {
  sf <- if (inherits(model, "opt_standard_form")) model else
    opt_standard_form(model)
  integer <- any(sf$type %in% c("integer", "binary"))
  rk <- sf$row_kind
  has_nl <- sf$objective_kind == "nonlinear" || any(rk == "nonlinear")
  has_qrow <- any(rk == "quadratic")
  ctypes <- cone_types(sf)
  conic <- length(ctypes) > 0L && !has_nl && !has_qrow
  base <- if (has_nl) "NLP" else if (has_qrow) "QCQP" else
    if (conic) (if (identical(ctypes, "soc")) "SOCP" else "CONIC") else
    if (sf$objective_kind == "quadratic") "QP" else "LP"
  # indicator and SOS constraints are discrete even without integer columns
  integer <- integer || length(sf$indicators) > 0L || length(sf$sos) > 0L
  cls <- if (integer) paste0("MI", base) else base
  evidence <- sprintf("syntactic: objective is %s; constraint rows: %s",
                      sf$objective_kind,
                      if (length(rk)) paste(names(table(rk)), table(rk),
                                            sep = " = ", collapse = ", ")
                      else "none")
  if (integer) {
    evidence <- c(evidence, sprintf("syntactic: %d integer or binary variables",
                                    sum(sf$type %in% c("integer", "binary"))))
  }
  if (length(ctypes)) {
    evidence <- c(evidence, sprintf(
      "structural: %d cone constraint(s) of type %s; cones are convex sets",
      length(sf$cones), paste(ctypes, collapse = ", ")))
  }
  if (length(sf$indicators) || length(sf$sos)) {
    evidence <- c(evidence, sprintf(
      "syntactic: %d indicator and %d SOS constraint(s), which are discrete",
      length(sf$indicators), length(sf$sos)))
  }
  convex <- "convex"
  obj_curv <- NULL
  if (has_nl) {
    convex <- "unknown"
    evidence <- c(evidence, paste(
      "unknown: convexity of general nonlinear expressions is not analyzed;",
      "solutions from local solvers are local optima"))
  } else {
    if (sf$objective_kind == "quadratic") {
      Qo <- if (sf$sense == "max") -sf$Q else sf$Q
      obj_curv <- psd_check(Qo, tol)
      what <- if (sf$sense == "min") "convex" else "concave"
      if (isTRUE(obj_curv$psd)) {
        evidence <- c(evidence, sprintf("%s: objective is %s (%s)",
                                        obj_curv$method, what, obj_curv$detail))
      } else if (isFALSE(obj_curv$psd)) {
        convex <- "nonconvex"
        evidence <- c(evidence, sprintf("%s: objective is not %s (%s)",
                                        obj_curv$method, what, obj_curv$detail))
      } else {
        convex <- "unknown"
        evidence <- c(evidence, sprintf("%s: %s", obj_curv$method,
                                        obj_curv$detail))
      }
    } else {
      evidence <- c(evidence, "structural: linear objective")
    }
    for (r in names(sf$row_Q)) {
      ri <- as.integer(r)
      Qr <- sf$row_Q[[r]]
      fin_up <- is.finite(sf$row_upper[ri])
      fin_lo <- is.finite(sf$row_lower[ri])
      lab <- sf$row_names[ri]
      if (fin_up && fin_lo) {
        convex <- "nonconvex"
        evidence <- c(evidence, sprintf(
          "syntactic: quadratic row %s is bounded on both sides (or an equality), which is nonconvex",
          lab))
        next
      }
      chk <- psd_check(if (fin_up) Qr else -Qr, tol)
      if (isFALSE(chk$psd)) {
        convex <- "nonconvex"
        evidence <- c(evidence, sprintf("%s: row %s defines a nonconvex set (%s)",
                                        chk$method, lab, chk$detail))
      } else if (is.na(chk$psd)) {
        if (convex == "convex") convex <- "unknown"
        evidence <- c(evidence, sprintf("%s: row %s, %s", chk$method, lab,
                                        chk$detail))
      } else {
        evidence <- c(evidence, sprintf("%s: row %s defines a convex set (%s)",
                                        chk$method, lab, chk$detail))
      }
    }
    if (!has_qrow && sf$objective_kind != "quadratic" && !length(ctypes)) {
      evidence <- c(evidence, "structural: all rows linear, feasible set is a polyhedron")
    }
  }
  if (integer && convex == "convex") {
    evidence <- c(evidence,
                  "the continuous relaxation is convex; the integer problem is not")
  }
  structure(list(class = cls, integer = integer, special = special_kinds(sf),
                 cone_types = ctypes,
                 objective = sf$objective_kind,
                 constraints = if (length(rk)) as.list(table(rk)) else list(),
                 convexity = convex, objective_curvature = obj_curv,
                 evidence = evidence), class = "opt_classification")
}

#' @export
print.opt_classification <- function(x, ...) {
  cat(sprintf("<opt_classification: %s, %s%s>\n", x$class, x$convexity,
              if (x$integer) " relaxation" else ""))
  cat(paste0("  - ", x$evidence), sep = "\n")
  invisible(x)
}
