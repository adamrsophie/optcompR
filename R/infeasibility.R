#' Explain why a model is infeasible
#'
#' Solves an elastic relaxation of the model: every constraint row gets
#' nonnegative slack variables that let it move toward feasibility, and with
#' `bounds = TRUE` the bounds of continuous and integer variables do too. The
#' relaxation minimizes the weighted total (L1) violation, so its solution
#' names rows and bounds that can move, and by how much, for the model to
#' become feasible at the least weighted total change. This is not the
#' smallest number of rows to change, nor an irreducible infeasible subset:
#' two copies of `x >= 1` against one `100 * x <= 0` relax the two copies
#' (total change 2) although relaxing the single upper row would suffice
#' (change 100). Special constraints (indicators and SOS), binary bounds and
#' integrality are kept as they are.
#'
#' The model itself is not changed; the relaxation is built on a copy (see
#' [opt_copy()]). For nonlinear models the relaxation is solved by a local
#' solver unless the model is proven convex, so the explanation is then a
#' local one.
#'
#' @param model An [opt_model()].
#' @param bounds Also relax variable bounds.
#' @param weights Optional named numeric vector of weights by constraint
#'   block name (and `".bounds"` for variable bounds); the default weight is
#'   1. A larger weight makes moving that block less attractive.
#' @param solver,... Passed to [opt_solve()] for the relaxation.
#' @param tol Violations at or below `tol` are reported as zero.
#' @return An `opt_infeasibility` list with `feasible` (whether the model is
#'   feasible, that is, no unweighted violation above `tol` is needed; the
#'   weights do not change this verdict), `total` (the weighted sum of
#'   violations), `violation` (their unweighted sum), `violations`, a data
#'   frame with columns `kind` (`"row"`
#'   or `"bound"`), `name`, `block`, `side` (`"lower"` or `"upper"`), `bound`
#'   (the original bound), `value` (the row or variable value at the
#'   relaxed solution) and `violation`, and `result`, the [opt_solve()]
#'   result of the relaxation.
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x", dim = 2, lower = 0)
#' opt_constraint(m, x[1] + x[2] <= 1, name = "capacity")
#' opt_constraint(m, x[1] >= 2, name = "demand_1")
#' opt_constraint(m, x[2] >= 0.5, name = "demand_2")
#' opt_explain_infeasibility(m)
#' @export
opt_explain_infeasibility <- function(model, bounds = TRUE, weights = NULL,
                                      solver = "auto", ..., tol = 1e-7) {
  check_model(model)
  if (!length(model$cons) && !bounds) {
    stop("The model has no constraint rows to relax.", call. = FALSE)
  }
  w_of <- function(key) {
    w <- if (!is.null(weights) && key %in% names(weights)) weights[[key]] else 1
    if (!is.numeric(w) || length(w) != 1L || !is.finite(w) || w <= 0) {
      stop("Weights must be positive finite numbers.", call. = FALSE)
    }
    w
  }
  em <- opt_copy(model)
  fresh <- function(base) {
    nm <- base
    k <- 1L
    while (!is.na(index_get(em$block_index, nm))) {
      k <- k + 1L
      nm <- paste0(base, "_", k)
    }
    nm
  }
  parts <- list()
  penalty <- 0
  for (k in seq_along(em$cons)) {
    blk <- em$cons[[k]]
    n <- blk$expr$len
    lo_ok <- is.finite(blk$lower)
    hi_ok <- is.finite(blk$upper)
    # up raises the row toward its lower bound, down lowers it to its upper
    up <- opt_variable(em, fresh(paste0("elastic_up_", blk$name)), dim = n,
                       lower = 0, upper = ifelse(lo_ok, Inf, 0))
    down <- opt_variable(em, fresh(paste0("elastic_down_", blk$name)),
                         dim = n, lower = 0, upper = ifelse(hi_ok, Inf, 0))
    blk$expr <- expr_add(blk$expr, expr_add(up, down, sign = -1))
    em$cons[[k]] <- blk
    w <- w_of(blk$name)
    penalty <- expr_add(penalty, expr_scale(expr_add(sum(up), sum(down)), w))
    parts[[length(parts) + 1L]] <- list(
      kind = "row", block = blk$name,
      names = element_labels(blk$name, n, model$cons[[k]]$expr$dim),
      expr = model$cons[[k]]$expr,
      lower = blk$lower, upper = blk$upper, up = up, down = down)
  }
  if (bounds) {
    # variables of the original model only, not the elastic ones
    vs <- which(model$sym_kind == "var" & model$type != "binary" &
                  (is.finite(model$lower) | is.finite(model$upper)))
    if (length(vs)) {
      lo <- em$lower[vs]
      hi <- em$upper[vs]
      em$lower[vs] <- -Inf
      em$upper[vs] <- Inf
      xv <- symbol_expr(em, vs, NULL)
      up <- opt_variable(em, fresh("elastic_up_bounds"), dim = length(vs),
                         lower = 0, upper = ifelse(is.finite(lo), Inf, 0))
      down <- opt_variable(em, fresh("elastic_down_bounds"), dim = length(vs),
                           lower = 0, upper = ifelse(is.finite(hi), Inf, 0))
      opt_constraint(em, expr_add(xv, expr_add(up, down, sign = -1)), lower = lo,
                     upper = hi, name = fresh("elastic_bounds"))
      w <- w_of(".bounds")
      penalty <- expr_add(penalty, expr_scale(expr_add(sum(up), sum(down)), w))
      parts[[length(parts) + 1L]] <- list(
        kind = "bound", block = NA_character_, names = em$sym_name[vs],
        expr = symbol_expr(model, vs, NULL), lower = lo, upper = hi,
        up = up, down = down)
    }
  }
  if (!is_expr(penalty)) {
    stop("The model has nothing to relax.", call. = FALSE)
  }
  opt_objective(em, penalty, sense = "min")
  touch_structure(em)
  r <- opt_solve(em, solver = solver, ...)
  if (!opt_status(r) %in% c("optimal", "locally_optimal", "almost_optimal")) {
    stop(sprintf(paste(
      "The elastic relaxation could not be solved (status '%s'). Integrality,",
      "binary bounds or special constraints alone may make the model",
      "infeasible."), opt_status(r)), call. = FALSE)
  }
  rows <- lapply(parts, function(p) {
    u <- as.numeric(opt_value(r, p$up))
    d <- as.numeric(opt_value(r, p$down))
    val <- eval_expr(p$expr, r$z)
    # report polynomial rows as written: constants moved to the bound
    k0 <- if (p$expr$type == "poly") {
      rep_len(compact_poly(p$expr)$const, length(val))
    } else 0
    val <- val - k0
    p$lower <- p$lower - k0
    p$upper <- p$upper - k0
    side <- ifelse(u > d, "lower", "upper")
    data.frame(kind = p$kind, name = p$names,
               block = if (is.na(p$block)) NA_character_ else p$block,
               side = side, bound = ifelse(u > d, p$lower, p$upper),
               value = val, violation = pmax(u, d), stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  tab$violation[tab$violation <= tol] <- 0
  tab <- tab[tab$violation > 0, , drop = FALSE]
  tab <- tab[order(-tab$violation), , drop = FALSE]
  rownames(tab) <- NULL
  total <- opt_objective_value(r)
  structure(list(feasible = nrow(tab) == 0L, total = total,
                 violation = sum(tab$violation), violations = tab,
                 result = r), class = "opt_infeasibility")
}

#' @export
print.opt_infeasibility <- function(x, ...) {
  if (x$feasible) {
    cat("<opt_infeasibility: the model is feasible; no constraint must move>\n")
    return(invisible(x))
  }
  v <- x$violations
  cat(sprintf(paste("<opt_infeasibility: %d row(s) and %d bound(s) must move,",
                    "total violation %s>\n"),
              sum(v$kind == "row"), sum(v$kind == "bound"), fmt_num(x$total)))
  n <- min(nrow(v), 20L)
  for (i in seq_len(n)) {
    cat(sprintf("  %s %s: %s bound %s, relaxed value %s (violation %s)\n", v$kind[i],
                v$name[i], v$side[i], fmt_num(v$bound[i]),
                fmt_num(v$value[i]), fmt_num(v$violation[i])))
  }
  if (nrow(v) > n) cat(sprintf("  ... %d more\n", nrow(v) - n))
  invisible(x)
}
