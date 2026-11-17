# Mixed-integer nonlinear programming by NLP-based branch and bound.
#
# Each node solves the continuous relaxation with NLopt over tightened
# bounds on the integer variables. NLopt stops on step and function
# tolerances, so its objective is a primal value, not a bound. When the
# relaxation is proven convex, each node gets a certified lower bound from
# an outer approximation instead: the objective linearized at the
# relaxation point (a supporting hyperplane of a convex function) minimized
# by HiGHS over the node's box and the linear rows. Only these bounds prune
# nodes and enter the reported bound, and a node closes only when its bound
# meets its best point; an integral point with a weak bound is branched
# further. When the relaxation is not known to be convex there are no
# bounds: the search is a heuristic and the result is locally optimal.

nlpbb_solve_sf <- function(sf, cls, control) {
  s <- if (sf$sense == "max") -1 else 1
  int <- which(sf$type %in% c("integer", "binary"))
  convex <- cls$convexity == "convex"
  relax <- cls
  relax$class <- sub("^MI", "", cls$class)
  relax$integer <- FALSE
  max_nodes <- control$options$max_nodes %||% 10000L
  int_tol <- control$options$integer_tolerance %||% 1e-6
  gap_tol <- control$mip_gap %||% 1e-9
  sub_opts <- control$options
  sub_opts$max_nodes <- NULL
  sub_opts$integer_tolerance <- NULL
  t_end <- proc.time()[["elapsed"]] + control$time_limit
  lower0 <- sf$lower
  upper0 <- sf$upper
  lower0[int] <- ceiling(lower0[int] - int_tol)
  upper0[int] <- floor(upper0[int] + int_tol)
  if (any(lower0 > upper0)) {
    return(nlpbb_result("infeasible", NULL, Inf, -Inf, 0L, convex, s, sf))
  }
  # open nodes: a list with a parallel vector of their bounds, the first
  # nopen entries in use; a node is removed by moving the last one into
  # its place, so nothing is copied or scanned in R code
  # ponytail: picking the best node is a which.min in C, linear in the open
  # nodes; a binary heap pays off only with very many cheap relaxations
  open <- list(list(lower = lower0, upper = upper0, bound = -Inf,
                    start = sf$start))
  ob <- -Inf
  nopen <- 1L
  incumbent <- Inf
  best <- NULL
  nodes <- 0L
  limit <- NULL
  # nodes whose relaxation stopped without a verdict (limits, numerical
  # trouble): their subtrees are lost, so nothing global can be claimed
  lost <- 0L
  # smallest relaxation bound of the subtrees pruned against the incumbent:
  # a subtree pruned only because of the requested gap keeps its bound, so
  # the final bound stays valid
  pruned <- Inf
  # smallest certified bound of the nodes closed as solved
  closed <- Inf
  while (nopen) {
    if (nodes >= max_nodes) { limit <- "node_limit"; break }
    left <- t_end - proc.time()[["elapsed"]]
    if (left <= 0) { limit <- "time_limit"; break }
    # best first: the open node with the smallest relaxation bound
    k <- which.min(ob[seq_len(nopen)])
    node <- open[[k]]
    open[k] <- open[nopen]
    ob[k] <- ob[nopen]
    open[nopen] <- list(NULL)
    nopen <- nopen - 1L
    if (convex && is.finite(incumbent) &&
        node$bound >= incumbent - gap_tol * max(1, abs(incumbent))) {
      pruned <- min(pruned, node$bound)
      next
    }
    nodes <- nodes + 1L
    sub <- sf
    sub$lower <- node$lower
    sub$upper <- node$upper
    sub$start <- node$start
    ctrl <- control
    ctrl$time_limit <- left
    ctrl$options <- sub_opts
    res <- tryCatch(nloptr_solve_sf(sub, relax, ctrl),
                    error = function(e) list(termination = "error",
                                             primal_status = "no_solution"))
    if (res$primal_status != "feasible_point" ||
        !res$termination %in% c("optimal", "locally_optimal")) {
      # a local NLP solver has no infeasibility certificate: a node that
      # looks infeasible may still hold solutions, so it counts as lost
      # unless the relaxation is convex and the solver says infeasible
      # "infeasible" comes only from the LP certificate in the NLopt
      # adapter; any other failure leaves the subtree unexplored
      if (!identical(res$termination, "infeasible")) {
        lost <- lost + 1L
      }
      next
    }
    val <- s * res$objective
    x <- res$x
    # certified lower bound of the node, never below its parent's
    lb <- -Inf
    if (convex) {
      oa <- oa_lower_bound(sub, x, s)
      lb <- max(node$bound, oa$bound)
      if (lb == Inf) next # box and linear rows infeasible
      # with linear rows only, the best cutting-plane point is a feasible
      # relaxation point; use it when NLopt stopped short of it
      if (all(sub$row_kind == "linear") && oa$f < val) {
        x <- oa$y
        val <- oa$f
      }
    }
    if (convex && is.finite(incumbent) &&
        lb >= incumbent - gap_tol * max(1, abs(incumbent))) {
      pruned <- min(pruned, lb)
      next
    }
    frac <- int[abs(x[int] - round(x[int])) > int_tol]
    if (!length(frac)) {
      if (val < incumbent) {
        # rounding within the tolerance can break other constraints, so
        # solve again with the integers fixed at their rounded values
        fix <- sub
        fix$lower[int] <- fix$upper[int] <- round(x[int])
        fix$start <- x
        fix$start[int] <- round(x[int])
        rf <- tryCatch(nloptr_solve_sf(fix, relax, ctrl),
                       error = function(e) NULL)
        if (!is.null(rf) && rf$primal_status == "feasible_point" &&
            rf$termination %in% c("optimal", "locally_optimal")) {
          fval <- s * rf$objective
          if (fval < incumbent) {
            incumbent <- fval
            best <- rf$x
            best[int] <- round(best[int])
          }
        } else {
          lost <- lost + 1L
        }
      }
      # the node is solved when its certified bound reaches its point;
      # otherwise a better integer point may remain in it, so it is split
      # at the middle of its widest integer range
      if (lb >= val - max(gap_tol, 1e-7) * max(1, abs(val))) {
        closed <- min(closed, lb)
        next
      }
      wide <- int[node$upper[int] > node$lower[int]]
      if (!length(wide)) {
        # all integers fixed and the bound still open: nothing to branch on
        lost <- lost + 1L
        next
      }
      j <- wide[which.max(node$upper[wide] - node$lower[wide])]
      cut <- floor((node$lower[j] + node$upper[j]) / 2)
    } else {
      j <- frac[which.max(abs(x[frac] - round(x[frac])))]
      cut <- floor(x[j])
    }
    dn <- node
    dn$upper[j] <- cut
    dn$bound <- lb
    dn$start <- pmin(pmax(x, dn$lower), dn$upper)
    up <- node
    up$lower[j] <- cut + 1
    up$bound <- lb
    up$start <- pmin(pmax(x, up$lower), up$upper)
    open[nopen + 1:2] <- list(dn, up)
    ob[nopen + 1:2] <- val
    nopen <- nopen + 2L
  }
  exact <- convex && lost == 0L
  bound <- if (!exact) -Inf else
    min(c(incumbent, pruned, closed, ob[seq_len(nopen)]))
  # with no incumbent, every node was certified infeasible or lost
  term <- if (!is.null(limit)) limit else if (is.null(best))
    (if (lost) "unknown" else "infeasible") else
    (if (exact) "optimal" else "locally_optimal")
  out <- nlpbb_result(term, best, incumbent, bound, nodes, convex, s, sf,
                      local = !exact)
  if (lost) {
    out$message <- sprintf(paste(
      "%s; %d node relaxation(s) stopped without a verdict, so the search",
      "is incomplete"), out$message, lost)
  }
  out
}

# Certified lower bound of the minimization form s * f of a node whose
# relaxation is convex. A convex f satisfies f(y) >= f(x_k) + g_k'(y - x_k)
# at every point x_k, and so does any convex combination h of these
# supporting hyperplanes. Kelley's method collects hyperplanes at the
# relaxation point and at successive LP solutions over the node's box and
# its linear rows. A nonlinear row side that composition rules prove convex
# on the node's box (g <= u with g convex, g >= l with g concave) adds its
# own supporting hyperplane at each point, an outer approximation of its
# set; other nonlinear rows are dropped, which only weakens the bound. The LP optimum itself is not used as the bound, since an LP
# solver stops within its tolerances; instead the LP multipliers give
#   min h(y) over l <= Ay <= u and the box
#     >= min over the box of (h(y) - lambda'Ay) + sum_i min over [l_i, u_i]
#        of lambda_i w_i,
# which holds for any lambda and is evaluated in closed form (weak duality
# with the box kept explicit). Returns the bound (-Inf when nothing
# finite is certified, Inf when the box and linear rows are infeasible) and
# the point with the smallest objective seen.
oa_lower_bound <- function(sf, x, s, iterations = 20L) {
  vs <- sf$var_sym
  n0 <- length(vs)
  n <- length(sf$c)
  ns <- sf$n_sym
  fg <- function(y) {
    z <- numeric(ns)
    z[vs] <- y[seq_len(n0)]
    jj <- tryCatch(jac_expr(sf$objective_lowered, z, ns),
                   error = function(e) NULL)
    if (is.null(jj)) return(NULL)
    out <- list(f = s * jj$val,
                g = c(s * as.vector(jj$J[, vs, drop = FALSE]), numeric(n - n0)))
    if (all(is.finite(c(out$f, out$g)))) out else NULL
  }
  lin <- which(sf$row_kind == "linear")
  A0 <- methods::as(sf$A[lin, , drop = FALSE], "CsparseMatrix")
  rl0 <- sf$row_lower[lin]
  ru0 <- sf$row_upper[lin]
  lo <- sf$lower
  up <- sf$upper
  # nonlinear row elements with a side proven convex on this box
  slo <- rep(-Inf, ns)
  shi <- rep(Inf, ns)
  slo[vs] <- lo[seq_len(n0)]
  shi[vs] <- up[seq_len(n0)]
  rowcut <- list()
  off <- 0L
  for (e in sf$cons_lowered) {
    if (is.null(e)) next
    rows <- off + seq_len(e$len)
    off <- off + e$len
    nonlin <- sf$row_kind[rows] != "linear"
    if (!any(nonlin)) next
    cv <- tryCatch(curvature(e, slo, shi)$curv, error = function(err) NULL)
    if (is.null(cv)) next
    hi_ok <- nonlin & is.finite(sf$expr_upper[rows]) &
      cv %in% c(CURV_CONST, CURV_AFFINE, CURV_CONVEX)
    lo_ok <- nonlin & is.finite(sf$expr_lower[rows]) &
      cv %in% c(CURV_CONST, CURV_AFFINE, CURV_CONCAVE)
    if (any(hi_ok | lo_ok)) {
      rowcut[[length(rowcut) + 1L]] <- list(e = e, hi = hi_ok, lo = lo_ok,
        u = sf$expr_upper[rows], l = sf$expr_lower[rows])
    }
  }
  # supporting hyperplanes of those rows at y: g(w) >= g(y) + J (w - y),
  # so g(w) <= u implies J w <= u - g(y) + J y (and >= for a concave side)
  row_cuts <- function(y) {
    z <- numeric(ns)
    z[vs] <- y[seq_len(n0)]
    C <- NULL
    l <- numeric()
    u <- numeric()
    for (rc in rowcut) {
      jj <- tryCatch(jac_expr(rc$e, z, ns), error = function(err) NULL)
      if (is.null(jj)) next
      J <- as.matrix(jj$J[, vs, drop = FALSE])
      if (n > n0) J <- cbind(J, matrix(0, nrow(J), n - n0))
      for (j in which(rc$hi | rc$lo)) {
        cst <- jj$val[j] - sum(J[j, ] * y)
        if (!all(is.finite(c(J[j, ], cst)))) next
        C <- rbind(C, J[j, ])
        l <- c(l, if (rc$lo[j]) rc$l[j] - cst else -Inf)
        u <- c(u, if (rc$hi[j]) rc$u[j] - cst else Inf)
      }
    }
    list(C = C, l = l, u = u)
  }
  A1 <- A0
  rl <- rl0
  ru <- ru0
  # min over [a, b] of r * y, exactly; -Inf when unbounded
  box_min <- function(r, a, b) {
    v <- ifelse(r > 0, r * a, ifelse(r < 0, r * b, 0))
    v[is.nan(v)] <- 0
    sum(v)
  }
  safe_bound <- function(G, b0, mu, lam) {
    mu <- pmax(mu, 0)
    if (!(sum(mu) > 0)) return(-Inf)
    mu <- mu / sum(mu)
    # h(y) = sum_k mu_k (f_k - g_k'x_k) + (sum_k mu_k g_k)'y
    hc <- as.vector(crossprod(G, mu))
    h0 <- sum(mu * b0)
    lam[(lam > 0 & !is.finite(rl)) | (lam < 0 & !is.finite(ru))] <- 0
    r <- hc - if (length(lam)) as.vector(Matrix::crossprod(A1, lam)) else 0
    wl <- ifelse(lam > 0, lam * rl, ifelse(lam < 0, lam * ru, 0))
    val <- h0 + box_min(r, lo, up) + sum(wl)
    # floating point error of the evaluation itself
    mag <- abs(h0) + sum(abs(r) * pmax(abs(lo), abs(up)), na.rm = TRUE) +
      sum(abs(wl))
    if (!is.finite(val)) return(-Inf)
    val - 1e-12 * max(1, mag)
  }
  G <- NULL
  b0 <- numeric()
  best <- list(f = Inf, y = x)
  lb <- -Inf
  pt <- x
  ctrl <- highs::highs_control(threads = .optcompR_state$highs_threads %||% 1L,
                               log_to_console = FALSE)
  for (it in seq_len(iterations)) {
    v <- fg(pt)
    if (is.null(v)) break
    if (v$f < best$f) best <- list(f = v$f, y = pt)
    # t >= f_k + g_k'(y - x_k), written g_k'y - t <= g_k'x_k - f_k
    G <- rbind(G, v$g)
    b0 <- c(b0, v$f - sum(v$g * pt))
    if (length(rowcut)) {
      rc <- row_cuts(pt)
      if (!is.null(rc$C)) {
        A1 <- rbind(A1, methods::as(unname(rc$C), "CsparseMatrix"))
        rl <- c(rl, rc$l)
        ru <- c(ru, rc$u)
      }
    }
    A <- rbind(cbind(A1, sparse_zero(length(rl), 1L)),
               methods::as(cbind(G, -1), "CsparseMatrix"))
    sol <- tryCatch({
      hm <- highs::highs_model(L = c(numeric(n), 1), lower = c(lo, -Inf),
                               upper = c(up, Inf), A = A,
                               lhs = c(rl, rep(-Inf, it)), rhs = c(ru, -b0))
      sv <- highs::highs_solver(hm, control = ctrl)
      sv$solve()
      list(status = sv$status_message(), sol = sv$solution())
    }, error = function(e) NULL)
    if (is.null(sol)) break
    if (identical(sol$status, "Infeasible") && it == 1L) {
      return(list(bound = Inf, y = best$y, f = best$f))
    }
    if (!identical(sol$status, "Optimal") || !isTRUE(sol$sol$dual_valid)) break
    rd <- sol$sol$row_dual
    nl <- length(rl)
    # cut rows are upper bounded: their multipliers are -row_dual
    lb <- max(lb, safe_bound(G, b0, -rd[nl + seq_len(it)], rd[seq_len(nl)]))
    pt <- pmin(pmax(sol$sol$col_value[seq_len(n)], lo), up)
    if (best$f - lb <= 1e-9 * max(1, abs(best$f))) break
  }
  list(bound = lb, y = best$y, f = best$f)
}

nlpbb_result <- function(term, best, incumbent, bound, nodes, convex, s, sf,
                         local = !convex) {
  has <- !is.null(best)
  gap <- if (has && is.finite(bound))
    abs(incumbent - bound) / max(1, abs(incumbent)) else NA_real_
  list(termination = term,
       primal_status = if (has) "feasible_point" else "no_solution",
       dual_status = "not_available",
       x = if (has) best else rep(NA_real_, length(sf$c)),
       objective = if (has) s * incumbent else NA_real_,
       row_dual = NULL, col_dual = NULL, iterations = c(nodes = nodes),
       gap = gap, bound = if (is.finite(bound)) s * bound else NA_real_,
       message = sprintf("NLP branch and bound, %d nodes, %s relaxation", nodes,
                         if (convex) "convex" else "not known to be convex"),
       local = local, warm_start = FALSE, extra = list(), raw = NULL)
}

register_nlpbb <- function() {
  register_solver(
    "nlpbb", "nloptr",
    "NLP branch and bound over NLopt relaxations (global for convex MINLP, heuristic otherwise)",
    classes = c("MINLP", "MIQCQP", "MIQP"), convex_only = character(),
    # a global search whose results mark themselves local when the
    # relaxation is not convex or the search is incomplete
    global = c("MINLP", "MIQCQP", "MIQP"),
    features = list(time_limit = TRUE, mip_gap = TRUE, warm_start = TRUE),
    solve = nlpbb_solve_sf, version = pkg_version("nloptr"), priority = 80L)
}

# Certify infeasibility when all rows are linear: an LP feasibility problem
# over the bounds, solved by HiGHS.
linear_infeasible <- function(sf) {
  if (!requireNamespace("highs", quietly = TRUE)) return(FALSE)
  if (any(sf$row_kind != "linear")) return(FALSE)
  if (any(sf$lower > sf$upper)) return(TRUE)
  # bounds alone are feasible exactly when every lower <= upper
  if (!nrow(sf$A)) return(FALSE)
  n <- length(sf$lower)
  A <- sf$A
  res <- tryCatch(highs::highs_solve(
    L = numeric(n), lower = sf$lower, upper = sf$upper,
    A = if (nrow(A)) A else NULL,
    lhs = if (nrow(A)) sf$row_lower else NULL,
    rhs = if (nrow(A)) sf$row_upper else NULL),
    error = function(e) NULL)
  !is.null(res) && identical(res$status_message, "Infeasible")
}
