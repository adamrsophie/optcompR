# Mixed-integer nonlinear programming by NLP-based branch and bound.
#
# Each node solves the continuous relaxation with NLopt over tightened
# bounds on the integer variables. When the relaxation is convex, node
# optima are lower bounds and pruning by bound is valid, so a completed
# search proves global optimality (up to the NLP tolerances). When the
# relaxation is not known to be convex, node optima are local, pruning can
# discard better solutions, and the result is reported as locally optimal:
# a heuristic, labeled as such.

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
  open <- list(list(lower = lower0, upper = upper0, bound = -Inf,
                    start = sf$start))
  incumbent <- Inf
  best <- NULL
  nodes <- 0L
  limit <- NULL
  # nodes whose relaxation stopped without a verdict (limits, numerical
  # trouble): their subtrees are lost, so nothing global can be claimed
  lost <- 0L
  while (length(open)) {
    if (nodes >= max_nodes) { limit <- "node_limit"; break }
    left <- t_end - proc.time()[["elapsed"]]
    if (left <= 0) { limit <- "time_limit"; break }
    # best first: the open node with the smallest relaxation bound
    k <- which.min(vapply(open, `[[`, 0, "bound"))
    node <- open[[k]]
    open[[k]] <- NULL
    if (convex && is.finite(incumbent) && node$bound >= incumbent - gap_tol * max(1, abs(incumbent))) next
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
    if (convex && is.finite(incumbent) && val >= incumbent - gap_tol * max(1, abs(incumbent))) next
    x <- res$x
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
      next
    }
    j <- frac[which.max(abs(x[frac] - round(x[frac])))]
    dn <- node
    dn$upper[j] <- floor(x[j])
    dn$bound <- val
    dn$start <- x
    up <- node
    up$lower[j] <- ceiling(x[j])
    up$bound <- val
    up$start <- x
    open <- c(open, list(dn, up))
  }
  exact <- convex && lost == 0L
  bound <- if (!exact) -Inf else if (length(open))
    min(c(incumbent, vapply(open, `[[`, 0, "bound"))) else incumbent
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
    features = list(time_limit = TRUE, mip_gap = TRUE, warm_start = TRUE),
    solve = nlpbb_solve_sf, version = pkg_version("nloptr"), priority = 80L)
}

# Certify infeasibility when all rows are linear: an LP feasibility problem
# over the bounds, solved by HiGHS.
linear_infeasible <- function(sf) {
  if (!requireNamespace("highs", quietly = TRUE)) return(FALSE)
  if (any(sf$row_kind != "linear")) return(FALSE)
  if (any(sf$lower > sf$upper)) return(TRUE)
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
