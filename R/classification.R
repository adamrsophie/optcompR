# Problem classification with recorded evidence.

# Curvature of a symmetric matrix: is it positive semidefinite?
# Returns list(psd = TRUE/FALSE/NA, method, detail).
psd_check <- function(Q, tol = 1e-9) {
  if (is.null(Q) || length(Q@x) == 0L || all(Q@x == 0)) {
    return(list(psd = TRUE, method = "structural", min_eigen = 0,
                detail = "zero matrix"))
  }
  Qs <- Matrix::forceSymmetric(Q)
  Tq <- methods::as(Qs, "TsparseMatrix")
  if (all(Tq@i == Tq@j)) {
    d <- Matrix::diag(Q)
    return(list(psd = all(d >= 0), method = "structural", min_eigen = min(d),
                detail = "diagonal matrix; signs of the diagonal decide exactly"))
  }
  # rows and columns without entries do not affect definiteness
  S <- sort(unique(c(Tq@i, Tq@j))) + 1L
  if (length(S) < nrow(Q)) {
    Qs <- Matrix::forceSymmetric(Qs[S, S, drop = FALSE])
  }
  n <- nrow(Qs)
  # a symmetric matrix with a nonnegative, weakly dominant diagonal is PSD
  # (Gershgorin), a cheap exact test before any factorization
  d <- Matrix::diag(Qs)
  off <- Matrix::rowSums(abs(Qs)) - abs(d)
  if (all(d >= off)) {
    return(list(psd = TRUE, method = "structural", min_eigen = NA_real_,
                detail = "diagonally dominant with a nonnegative diagonal (Gershgorin)"))
  }
  if (n <= 2000L) {
    ev <- eigen(as.matrix(Qs), symmetric = TRUE, only.values = TRUE)$values
    # relative to the largest eigenvalue, so tiny matrices are not
    # declared convex by an absolute tolerance
    scale <- max(abs(ev))
    lo <- min(ev)
    # eigenvalue rounding error is about n * eps * ||Q||; a negative
    # eigenvalue beyond it but inside `tol` cannot be told from real
    # negative curvature, so convexity is left undetermined
    round_tol <- 100 * n * .Machine$double.eps * scale
    status <- if (lo >= -round_tol) TRUE else if (lo >= -tol * scale) NA else
      FALSE
    detail <- sprintf("smallest eigenvalue %s (rounding level %s, tolerance %s)",
                      fmt_num(lo), fmt_num(round_tol), fmt_num(tol * scale))
    if (isTRUE(status) && lo < round_tol) {
      detail <- paste0(detail, "; zero up to rounding")
    }
    if (is.na(status)) {
      detail <- paste0(detail, "; slightly negative, so curvature is undetermined")
    }
    return(list(psd = status, method = "numerical eigenvalues",
                min_eigen = lo, detail = detail))
  }
  # a shift at rounding level only: a larger one would prove Q + shift I
  # positive definite while Q itself has negative curvature
  shift <- 100 * n * .Machine$double.eps * max(abs(Q@x))
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
#' that cannot be established automatically. Nonlinear expressions are
#' analyzed by composition rules (disciplined convex programming) on the box
#' of variable bounds: for example `exp()` of a convex expression is convex,
#' `log()` and `sqrt()` of a concave one are concave, and `sin()` is concave
#' where its argument lies in \eqn{[0, \pi]}. When the rules prove every
#' part, the model is `"convex"`; otherwise it is `"unknown"`, never guessed
#' `"nonconvex"`.
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
    dcp <- tryCatch(dcp_convexity(sf, tol), error = function(e) NULL)
    convex <- if (isTRUE(dcp$convex)) "convex" else "unknown"
    evidence <- c(evidence, dcp$evidence %||% paste(
      "unknown: convexity of these nonlinear expressions could not be",
      "analyzed; solutions from local solvers are local optima"))
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
    qg <- row_quad_groups(sf)
    rq <- sf$row_quad
    # evidence per row is collected in a preallocated vector, and each row's
    # matrix lives on its own support only
    rev_ <- character(length(qg))
    Ar <- NULL
    for (g in seq_along(qg)) {
      ri <- as.integer(names(qg)[g])
      qi <- qg[[g]]
      # a separable row (squares only) is decided by the signs of its
      # summed diagonal terms, without building its matrix
      sep <- all(rq$j[qi] == rq$k[qi])
      fin_up <- is.finite(sf$row_upper[ri])
      fin_lo <- is.finite(sf$row_lower[ri])
      lab <- sf$row_names[ri]
      if (!fin_up && !fin_lo) next # no bound, so no restriction
      if (fin_up && fin_lo) {
        # a side that the quadratic can never cross is redundant: for
        # 0 <= x^2 <= 1 only the upper side restricts, and the set is convex
        if (is.null(Ar)) Ar <- methods::as(sf$A, "RsparseMatrix")
        loc <- row_local(rq, qi, Ar, ri)
        two <- two_sided_quadratic(loc$Q, loc$a, sf$row_lower[ri],
                                   sf$row_upper[ri])
        if (two$verdict == "nonconvex") convex <- "nonconvex" else
          if (two$verdict == "unknown" && convex == "convex") convex <- "unknown"
        rev_[g] <- sprintf("%s: quadratic row %s is bounded on both sides; %s",
                           two$method, lab, two$detail)
        next
      }
      chk <- if (sep) {
        d <- rowsum(rq$v[qi], rq$j[qi])[, 1L] * (if (fin_up) 1 else -1)
        list(psd = all(d >= 0), method = "structural",
             detail = "diagonal matrix; signs of the diagonal decide exactly")
      } else {
        Ql <- row_local(rq, qi)$Q
        psd_check(if (fin_up) Ql else -Ql, tol)
      }
      if (isFALSE(chk$psd)) {
        convex <- "nonconvex"
        rev_[g] <- sprintf("%s: row %s defines a nonconvex set (%s)",
                           chk$method, lab, chk$detail)
      } else if (is.na(chk$psd)) {
        if (convex == "convex") convex <- "unknown"
        rev_[g] <- sprintf("%s: row %s, %s", chk$method, lab, chk$detail)
      } else {
        rev_[g] <- sprintf("%s: row %s defines a convex set (%s)",
                           chk$method, lab, chk$detail)
      }
    }
    evidence <- c(evidence, rev_[nzchar(rev_)])
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

# The quadratic and linear parts of quadratic row r on its own support: the
# columns named by either index of a quadratic term (a symmetric matrix keeps
# one triangle, so both are needed) and, when the row-compressed constraint
# matrix Ar is given, by a linear coefficient. Q is the |S| x |S| matrix of
# the one-half convention, a the linear coefficients on S.
row_local <- function(rq, qi, Ar = NULL, r = NULL) {
  j <- rq$j[qi]
  k <- rq$k[qi]
  v <- rq$v[qi]
  lj <- integer()
  lv <- numeric()
  if (!is.null(Ar)) {
    w <- seq.int(Ar@p[r] + 1L, length.out = Ar@p[r + 1L] - Ar@p[r])
    lj <- Ar@j[w] + 1L
    lv <- Ar@x[w]
    # a stored zero (a parameter at 0) is no support
    keep <- lv != 0
    lj <- lj[keep]
    lv <- lv[keep]
  }
  S <- sort(unique(c(j, k, lj)))
  Q <- half_q(match(j, S), match(k, S), v, length(S))
  a <- numeric(length(S))
  a[match(lj, S)] <- lv
  list(S = S, Q = Q, a = a)
}

# Convexity of the set lo <= a'x + 0.5 x'Qx <= up with both bounds finite,
# for Q and a on the row's own support. When Q is PSD (or NSD) and the
# quadratic attains its minimum (maximum) m, the lower (upper) side is
# redundant if lo <= m (up >= m), and the set is then the convex sublevel
# (superlevel) set. Otherwise, for Q PSD, nonzero and m < lo, a line through
# the minimizer along a curved direction meets the set in two separate
# intervals, so the set is nonconvex. Everything else is left unknown.
#
# Redundancy changes the topology of the set, so it is accepted only when
# lo <= m holds beyond the rounding error of computing m; a lower bound a
# little above the minimum (1e-5 <= x^2) is never treated as redundant,
# whatever feasibility tolerance a solver uses.
two_sided_quadratic <- function(Q, a, lo, up) {
  n <- length(a)
  if (n > 500L) {
    return(list(verdict = "unknown", method = "numerical",
                detail = "too many variables to decide redundancy of a side"))
  }
  Qs <- as.matrix(Q)
  e <- eigen(Qs, symmetric = TRUE)
  scale <- max(1, abs(e$values))
  eps <- 100 * n * .Machine$double.eps
  rt <- eps * scale
  for (sg in c(1, -1)) {
    ev <- sg * e$values
    if (min(ev) < -rt) next
    # extremum of sg * g: solve sg * (Q x + a) = 0 by the pseudo-inverse
    keep <- ev > rt
    V <- e$vectors[, keep, drop = FALSE]
    xs <- -V %*% ((crossprod(V, a)) / e$values[keep])
    # attained only if the gradient vanishes up to rounding
    attained <- max(abs(Qs %*% xs + a)) <= eps * (scale * max(1, abs(xs)) +
                                                    max(1, abs(a)))
    if (!attained) next
    m <- sum(a * xs) + 0.5 * sum(xs * (Qs %*% xs))
    # rounding error of m
    err <- eps * (abs(m) + sum(abs(a * xs)) + scale * sum(xs^2))
    if ((sg > 0 && lo <= m - err) || (sg < 0 && up >= m + err) ||
        (sg > 0 && lo == m) || (sg < 0 && up == m)) {
      return(list(verdict = "convex", method = "numerical", detail = sprintf(
        "the %s side is redundant (the quadratic's %s is %s), so the set is convex",
        if (sg > 0) "lower" else "upper", if (sg > 0) "minimum" else "maximum",
        fmt_num(m))))
    }
    if (any(keep) && ((sg > 0 && lo > m + err) || (sg < 0 && up < m - err))) {
      return(list(verdict = "nonconvex", method = "numerical", detail = sprintf(
        "the %s %s lies outside the bounds, so the set is not convex",
        if (sg > 0) "minimum" else "maximum", fmt_num(m))))
    }
  }
  list(verdict = "unknown", method = "numerical",
       detail = "neither side is shown redundant and nonconvexity is not proven")
}
