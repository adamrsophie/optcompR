# Exact lifting of polynomial nonlinear models to quadratic ones, for SCIP.
#
# A product of two factors whose degrees add up to more than two gets
# auxiliary variables w = factor for the quadratic factors, with the
# quadratic equality rows w - factor = 0. Integer powers use repeated
# squaring. The lifted model has the same feasible set and optimum in the
# original variables; only +, *, linear maps, division by constants and
# nonnegative integer powers can be lifted. Auxiliary variables get bounds
# from interval arithmetic on the variable bounds.

# Why a standard form cannot be lifted, or NULL when it can.
lift_refusal <- function(sf) {
  bad <- character()
  walk <- function(e) {
    if (e$type == "poly") return(invisible())
    op <- e$op
    ok <- op %in% c("+", "*", "linmap") ||
      (op == "/" && e$args[[2L]]$type == "poly" &&
         !length(e$args[[2L]]$lv) && !length(e$args[[2L]]$qv) &&
         all(e$args[[2L]]$const != 0)) ||
      (op == "^" && all(e$data$p >= 0 & e$data$p == round(e$data$p)))
    if (!ok) bad <<- c(bad, if (op == "fn") "external functions" else
      if (op == "^") "non-integer or negative powers" else
        if (op == "/") "division by expressions" else op)
    for (a in e$args) walk(a)
  }
  if (!is.null(sf$objective_lowered)) walk(sf$objective_lowered)
  for (e in sf$cons_lowered) if (!is.null(e)) walk(e)
  if (!length(bad)) return(NULL)
  sprintf(paste("scip handles polynomial nonlinear models only (lifted",
                "exactly to quadratic rows); this model uses %s"),
          paste(unique(bad), collapse = ", "))
}

lift_sf <- function(sf) {
  ns <- sf$n_sym
  n <- length(sf$c)
  st <- new.env(parent = emptyenv())
  st$naux <- 0L
  st$defs <- list()
  slo <- rep(-Inf, ns)
  shi <- rep(Inf, ns)
  slo[sf$var_sym] <- sf$lower
  shi[sf$var_sym] <- sf$upper
  st$slo <- slo
  st$shi <- shi
  # w = p for every element of the polynomial p
  aux_of <- function(p) {
    k <- p$len
    ids <- ns + st$naux + seq_len(k)
    iv <- curvature(p, st$slo, st$shi)
    st$naux <- st$naux + k
    st$slo <- c(st$slo, iv$lo)
    st$shi <- c(st$shi, iv$hi)
    st$defs[[length(st$defs) + 1L]] <- list(ids = ids, poly = p)
    new_poly(NULL, k, NULL, const = numeric(k), li = seq_len(k), lj = ids,
             lv = rep.int(1, k))
  }
  affine <- function(p) if (poly_degree(p) <= 1L) p else aux_of(p)
  mul <- function(a, b) {
    a <- compact_poly(a)
    b <- compact_poly(b)
    if (poly_degree(a) + poly_degree(b) > 2L) {
      same <- identical(a, b)
      a <- affine(a)
      b <- if (same) a else affine(b)
    }
    expr_mul(a, b)
  }
  ipow <- function(a, p) {
    out <- new_poly(NULL, a$len, NULL, const = rep(1, a$len))
    base <- a
    while (p > 0) {
      if (p %% 2 == 1) out <- mul(out, base)
      p <- p %/% 2
      if (p > 0) base <- mul(base, base)
    }
    out
  }
  lift <- function(e) {
    if (e$type == "poly") return(e)
    a <- lapply(e$args, lift)
    switch(e$op,
      "+" = expr_add(a[[1L]], a[[2L]]),
      "*" = mul(a[[1L]], a[[2L]]),
      linmap = expr_linmap(e$data$M, a[[1L]]),
      "/" = expr_scale(a[[1L]], 1 / rep_len(a[[2L]]$const, e$len)),
      "^" = {
        p <- rep_len(e$data$p, e$len)
        b <- if (a[[1L]]$len == e$len) a[[1L]] else
          expr_subset(a[[1L]], rep.int(1L, e$len))
        if (length(unique(p)) == 1L) ipow(b, p[1L]) else {
          parts <- lapply(seq_len(e$len), function(i)
            ipow(expr_subset(b, i), p[i]))
          Reduce(function(acc, i) expr_add(acc, expr_linmap(
            Matrix::sparseMatrix(i = i, j = 1L, x = 1, dims = c(e$len, 1L)),
            parts[[i]])), seq_len(e$len), expr_scale(b, 0))
        }
      },
      stop(sprintf("Cannot lift operator '%s'.", e$op), call. = FALSE))
  }
  obj <- if (identical(sf$objective_kind, "nonlinear"))
    compact_poly(lift(sf$objective_lowered)) else NULL
  cons <- lapply(sf$cons_lowered, function(e)
    if (is.null(e)) NULL else compact_poly(lift(e)))
  na <- st$naux
  N <- n + na
  col <- function(s) ifelse(s <= ns, sf$col_of_sym[pmin(s, ns)], n + s - ns)
  # terms of element i of polynomial p, in columns
  terms <- function(p, i) {
    w <- which(p$li == i)
    q <- which(p$qi == i)
    list(k = rep_len(p$const, p$len)[i], j = col(p$lj[w]), v = p$lv[w],
         Q = if (length(q)) half_q(col(p$qj[q]), col(p$qk[q]), p$qv[q], N)
         else NULL)
  }
  pad <- function(M) {
    T <- methods::as(M, "TsparseMatrix")
    Matrix::sparseMatrix(i = T@i + 1L, j = T@j + 1L, x = T@x,
                         dims = c(nrow(M), N))
  }
  out <- sf
  A <- methods::as(pad(sf$A), "TsparseMatrix")
  ai <- A@i + 1L
  aj <- A@j + 1L
  ax <- A@x
  # symmetric matrices padded to N columns, both triangles kept
  pad_sym <- function(Q) {
    T <- methods::as(methods::as(Q, "generalMatrix"), "TsparseMatrix")
    Matrix::sparseMatrix(i = T@i + 1L, j = T@j + 1L, x = T@x, dims = c(N, N))
  }
  row_Q <- lapply(sf$row_Q, pad_sym)
  lo <- sf$row_lower
  up <- sf$row_upper
  kinds <- sf$row_kind
  for (r in which(sf$row_kind == "nonlinear")) {
    t <- terms(cons[[sf$row_block[r]]], sf$row_element[r])
    keep <- ai != r
    ai <- c(ai[keep], rep.int(r, length(t$j)))
    aj <- c(aj[keep], t$j)
    ax <- c(ax[keep], t$v)
    if (!is.null(t$Q)) row_Q[[as.character(r)]] <- t$Q
    lo[r] <- sf$expr_lower[r] - t$k
    up[r] <- sf$expr_upper[r] - t$k
    kinds[r] <- if (is.null(t$Q)) "linear" else "quadratic"
  }
  m0 <- nrow(sf$A)
  r <- m0
  for (d in st$defs) {
    for (i in seq_along(d$ids)) {
      r <- r + 1L
      t <- terms(d$poly, i)
      ai <- c(ai, rep.int(r, length(t$j) + 1L))
      aj <- c(aj, col(d$ids[i]), t$j)
      ax <- c(ax, 1, -t$v)
      if (!is.null(t$Q)) row_Q[[as.character(r)]] <- -t$Q
      lo <- c(lo, t$k)
      up <- c(up, t$k)
      kinds <- c(kinds, if (is.null(t$Q)) "linear" else "quadratic")
    }
  }
  out$A <- Matrix::sparseMatrix(i = ai, j = aj, x = ax, dims = c(r, N))
  out$row_Q <- row_Q
  out$row_lower <- lo
  out$row_upper <- up
  out$row_kind <- kinds
  out$row_names <- c(sf$row_names, sprintf("lift_%d", seq_len(r - m0)))
  if (!is.null(obj)) {
    t <- terms(obj, 1L)
    cc <- numeric(N)
    cc[t$j] <- cc[t$j] + t$v
    out$c <- cc
    out$Q <- t$Q
    out$offset <- sf$offset + t$k
    out$objective_kind <- if (is.null(t$Q)) "linear" else "quadratic"
  } else {
    out$c <- c(sf$c, numeric(na))
    if (!is.null(sf$Q)) out$Q <- pad_sym(sf$Q)
  }
  out$lower <- c(sf$lower, st$slo[ns + seq_len(na)])
  out$upper <- c(sf$upper, st$shi[ns + seq_len(na)])
  out$type <- c(sf$type, rep("continuous", na))
  out$start <- c(sf$start, rep(NA_real_, na))
  out$col_names <- c(sf$col_names, sprintf("lift_aux_%d", seq_len(na)))
  out$n_model <- n
  out
}
