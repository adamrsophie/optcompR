# Curvature of nonlinear expressions by composition rules (disciplined
# convex programming), with interval bounds of every subexpression from the
# variable bounds. A proof is sound: "convex" means each element is convex
# on the box of variable bounds; failing to prove it gives "unknown", never
# "nonconvex".

# curvature codes
CURV_CONST <- 0L
CURV_AFFINE <- 1L
CURV_CONVEX <- 2L
CURV_CONCAVE <- 3L
CURV_UNKNOWN <- 4L

curv_add <- function(a, b) {
  out <- pmax(a, b)
  out[a == CURV_CONST] <- b[a == CURV_CONST]
  out[b == CURV_CONST] <- a[b == CURV_CONST]
  both <- a != CURV_CONST & b != CURV_CONST
  same <- a == b | a == CURV_AFFINE | b == CURV_AFFINE
  out[both & !same] <- CURV_UNKNOWN
  out[both & same] <- pmax(a, b)[both & same]
  out
}

curv_scale <- function(cv, s) {
  out <- cv
  neg <- s < 0
  out[neg & cv == CURV_CONVEX] <- CURV_CONCAVE
  out[neg & cv == CURV_CONCAVE] <- CURV_CONVEX
  out[s == 0] <- CURV_CONST
  out
}

# interval products with 0 * Inf = 0
iv_mul <- function(alo, ahi, blo, bhi) {
  p <- cbind(alo * blo, alo * bhi, ahi * blo, ahi * bhi)
  p[is.nan(p)] <- 0
  list(lo = apply(p, 1L, min), hi = apply(p, 1L, max))
}

# Curvature and interval of each element of a lowered expression; slo and
# shi are bounds per symbol.
curvature <- function(e, slo, shi, tol = 1e-9) {
  n <- e$len
  rep_n <- function(v) rep_len(v, n)
  if (e$type == "poly") {
    lo <- e$const
    hi <- e$const
    if (length(e$lv)) {
      t <- iv_mul(e$lv, e$lv, slo[e$lj], shi[e$lj])
      lo <- lo + acc_sum(e$li, t$lo, n)
      hi <- hi + acc_sum(e$li, t$hi, n)
    }
    cv <- ifelse(acc_sum(e$li, abs(e$lv), n) > 0, CURV_AFFINE, CURV_CONST)
    if (length(e$qv)) {
      sq <- e$qj == e$qk
      t <- iv_mul(slo[e$qj], shi[e$qj], slo[e$qk], shi[e$qk])
      # a square is nonnegative
      t$lo[sq] <- ifelse(slo[e$qj[sq]] <= 0 & shi[e$qj[sq]] >= 0, 0,
                         pmin(slo[e$qj[sq]]^2, shi[e$qj[sq]]^2))
      t <- iv_mul(e$qv, e$qv, t$lo, t$hi)
      lo <- lo + acc_sum(e$qi, t$lo, n)
      hi <- hi + acc_sum(e$qi, t$hi, n)
      for (k in unique(e$qi)) {
        w <- e$qi == k
        sy <- sort(unique(c(e$qj[w], e$qk[w])))
        Q <- half_q(match(e$qj[w], sy), match(e$qk[w], sy), e$qv[w],
                    length(sy))
        up <- psd_check(Q, tol)$psd
        cv[k] <- if (isTRUE(up)) CURV_CONVEX else
          if (isTRUE(psd_check(-Q, tol)$psd)) CURV_CONCAVE else CURV_UNKNOWN
      }
    }
    return(list(curv = cv, lo = lo, hi = hi))
  }
  A <- lapply(e$args, curvature, slo = slo, shi = shi, tol = tol)
  if (e$op == "linmap") {
    M <- e$data$M
    a <- A[[1L]]
    Mp <- M
    Mp@x <- pmax(Mp@x, 0)
    Mn <- M
    Mn@x <- pmin(Mn@x, 0)
    fin <- function(v, s) {
      v[is.nan(v)] <- s
      v
    }
    lo <- fin(as.vector(Mp %*% a$lo + Mn %*% a$hi), -Inf)
    hi <- fin(as.vector(Mp %*% a$hi + Mn %*% a$lo), Inf)
    Tm <- methods::as(M, "TsparseMatrix")
    cv <- rep(CURV_CONST, n)
    for (t in seq_along(Tm@x)) {
      r <- Tm@i[t] + 1L
      cv[r] <- curv_add(cv[r], curv_scale(a$curv[Tm@j[t] + 1L], Tm@x[t]))
    }
    return(list(curv = cv, lo = lo, hi = hi))
  }
  if (e$op == "fn") {
    return(list(curv = rep(CURV_UNKNOWN, n), lo = rep(-Inf, n),
                hi = rep(Inf, n)))
  }
  if (length(A) == 2L) {
    a <- lapply(A[[1L]], rep_n)
    b <- lapply(A[[2L]], rep_n)
    a$curv <- as.integer(a$curv)
    b$curv <- as.integer(b$curv)
    if (e$op == "+") {
      return(list(curv = curv_add(a$curv, b$curv), lo = a$lo + b$lo,
                  hi = a$hi + b$hi))
    }
    if (e$op == "*") {
      r <- iv_mul(a$lo, a$hi, b$lo, b$hi)
      cv <- rep(CURV_UNKNOWN, n)
      ca <- a$curv == CURV_CONST
      cb <- b$curv == CURV_CONST
      cv[ca] <- curv_scale(b$curv[ca], a$lo[ca])
      cv[cb] <- curv_scale(a$curv[cb], b$lo[cb])
      cv[ca & cb] <- CURV_CONST
      return(list(curv = cv, lo = r$lo, hi = r$hi))
    }
    if (e$op == "/") {
      pos <- b$lo > 0
      neg <- b$hi < 0
      ilo <- ifelse(pos | neg, 1 / b$hi, -Inf)
      ihi <- ifelse(pos | neg, 1 / b$lo, Inf)
      r <- iv_mul(a$lo, a$hi, ilo, ihi)
      cv <- rep(CURV_UNKNOWN, n)
      cb <- b$curv == CURV_CONST
      cv[cb] <- curv_scale(a$curv[cb], 1 / b$lo[cb])
      # c / g: 1/g is convex and decreasing for g > 0, concave and
      # decreasing for g < 0
      ca <- a$curv == CURV_CONST & !cb
      inv <- atom_rule(b$curv, ifelse(pos, CURV_CONVEX, ifelse(neg, CURV_CONCAVE,
                                                               CURV_UNKNOWN)),
                       ifelse(pos | neg, -1L, 0L))
      cv[ca] <- curv_scale(inv[ca], a$lo[ca])
      return(list(curv = cv, lo = r$lo, hi = r$hi))
    }
    return(list(curv = rep(CURV_UNKNOWN, n), lo = rep(-Inf, n),
                hi = rep(Inf, n)))
  }
  g <- lapply(A[[1L]], rep_n)
  g$curv <- as.integer(g$curv)
  if (e$op == "^") {
    p <- rep_n(e$data$p)
    return(power_rule(g, p))
  }
  unary_rule(e$op, g)
}

# DCP composition: atom curvature f (convex/concave/unknown) and
# monotonicity mono (1 nondecreasing, -1 nonincreasing, 0 neither) on the
# range of g.
atom_rule <- function(gcurv, f, mono) {
  n <- length(gcurv)
  f <- rep_len(f, n)
  mono <- rep_len(mono, n)
  out <- rep(CURV_UNKNOWN, n)
  aff <- gcurv <= CURV_AFFINE
  vex <- f == CURV_CONVEX & (aff | (mono > 0 & gcurv == CURV_CONVEX) |
                               (mono < 0 & gcurv == CURV_CONCAVE))
  cav <- f == CURV_CONCAVE & (aff | (mono > 0 & gcurv == CURV_CONCAVE) |
                                (mono < 0 & gcurv == CURV_CONVEX))
  out[vex] <- CURV_CONVEX
  out[cav] <- CURV_CONCAVE
  out[gcurv == CURV_CONST] <- CURV_CONST
  out
}

unary_rule <- function(op, g) {
  lo <- g$lo
  hi <- g$hi
  n <- length(lo)
  pos <- lo >= 0
  negv <- hi <= 0
  U <- CURV_UNKNOWN
  r <- switch(op,
    exp = list(f = CURV_CONVEX, m = 1L, lo = exp(lo), hi = exp(hi)),
    expm1 = list(f = CURV_CONVEX, m = 1L, lo = expm1(lo), hi = expm1(hi)),
    log = list(f = ifelse(lo > 0, CURV_CONCAVE, U), m = 1L,
               lo = suppressWarnings(log(pmax(lo, 0))),
               hi = suppressWarnings(log(pmax(hi, 0)))),
    log1p = list(f = ifelse(lo > -1, CURV_CONCAVE, U), m = 1L,
                 lo = suppressWarnings(log1p(pmax(lo, -1))),
                 hi = suppressWarnings(log1p(pmax(hi, -1)))),
    sqrt = list(f = ifelse(lo >= 0, CURV_CONCAVE, U), m = 1L,
                lo = sqrt(pmax(lo, 0)), hi = sqrt(pmax(hi, 0))),
    abs = list(f = CURV_CONVEX, m = ifelse(pos, 1L, ifelse(negv, -1L, 0L)),
               lo = ifelse(pos, lo, ifelse(negv, -hi, 0)),
               hi = pmax(abs(lo), abs(hi))),
    cosh = list(f = CURV_CONVEX, m = ifelse(pos, 1L, ifelse(negv, -1L, 0L)),
                lo = ifelse(pos, cosh(lo), ifelse(negv, cosh(hi), 1)),
                hi = pmax(cosh(lo), cosh(hi))),
    sinh = list(f = ifelse(pos, CURV_CONVEX, ifelse(negv, CURV_CONCAVE, U)),
                m = 1L, lo = sinh(lo), hi = sinh(hi)),
    tanh = list(f = ifelse(pos, CURV_CONCAVE, ifelse(negv, CURV_CONVEX, U)),
                m = 1L, lo = tanh(lo), hi = tanh(hi)),
    atan = list(f = ifelse(pos, CURV_CONCAVE, ifelse(negv, CURV_CONVEX, U)),
                m = 1L, lo = atan(lo), hi = atan(hi)),
    asin = list(f = ifelse(pos & hi <= 1, CURV_CONVEX,
                           ifelse(negv & lo >= -1, CURV_CONCAVE, U)),
                m = 1L, lo = asin(pmax(lo, -1)), hi = asin(pmin(hi, 1))),
    acos = list(f = ifelse(pos & hi <= 1, CURV_CONCAVE,
                           ifelse(negv & lo >= -1, CURV_CONVEX, U)),
                m = -1L, lo = acos(pmin(hi, 1)), hi = acos(pmax(lo, -1))),
    sin = list(f = ifelse(pos & hi <= pi, CURV_CONCAVE,
                          ifelse(lo >= -pi & negv, CURV_CONVEX, U)),
               m = ifelse(lo >= -pi / 2 & hi <= pi / 2, 1L, 0L),
               lo = rep(-1, n), hi = rep(1, n)),
    cos = list(f = ifelse(lo >= -pi / 2 & hi <= pi / 2, CURV_CONCAVE,
                          ifelse(lo >= pi / 2 & hi <= 3 * pi / 2, CURV_CONVEX,
                                 U)),
               m = ifelse(pos & hi <= pi, -1L, ifelse(negv & lo >= -pi, 1L, 0L)),
               lo = rep(-1, n), hi = rep(1, n)),
    tan = list(f = ifelse(pos & hi < pi / 2, CURV_CONVEX,
                          ifelse(negv & lo > -pi / 2, CURV_CONCAVE, U)),
               m = ifelse(lo > -pi / 2 & hi < pi / 2, 1L, 0L),
               lo = ifelse(lo > -pi / 2 & hi < pi / 2, tan(lo), -Inf),
               hi = ifelse(lo > -pi / 2 & hi < pi / 2, tan(hi), Inf)),
    list(f = U, m = 0L, lo = rep(-Inf, n), hi = rep(Inf, n)))
  r$lo[is.nan(r$lo)] <- -Inf
  r$hi[is.nan(r$hi)] <- Inf
  list(curv = atom_rule(g$curv, r$f, r$m), lo = r$lo, hi = r$hi)
}

power_rule <- function(g, p) {
  lo <- g$lo
  hi <- g$hi
  n <- length(lo)
  pos <- lo >= 0
  negv <- hi <= 0
  U <- CURV_UNKNOWN
  even <- p > 0 & p == round(p) & p %% 2 == 0
  odd <- p > 1 & p == round(p) & p %% 2 == 1
  f <- ifelse(p == 1, CURV_AFFINE,
       ifelse(even, CURV_CONVEX,
       ifelse(p > 1, ifelse(pos, CURV_CONVEX,
                            ifelse(odd & negv, CURV_CONCAVE, U)),
       ifelse(p > 0, ifelse(pos, CURV_CONCAVE, U),
       ifelse(p < 0, ifelse(lo > 0, CURV_CONVEX, U), CURV_CONST)))))
  m <- ifelse(even, ifelse(pos, 1L, ifelse(negv, -1L, 0L)),
       ifelse(p > 0, 1L, ifelse(p < 0, -1L, 0L)))
  a <- suppressWarnings(lo^p)
  b <- suppressWarnings(hi^p)
  plo <- ifelse(even & !pos & !negv, 0, pmin(a, b))
  phi <- ifelse(even & !pos & !negv, pmax(a, b), pmax(a, b))
  plo[is.nan(plo)] <- -Inf
  phi[is.nan(phi)] <- Inf
  cv <- atom_rule(g$curv, f, m)
  cv[p == 1] <- g$curv[p == 1]
  cv[p == 0] <- CURV_CONST
  list(curv = cv, lo = plo, hi = phi)
}

# Prove convexity of a nonlinear standard form by composition rules.
# Returns list(convex = TRUE/FALSE, evidence = character()).
dcp_convexity <- function(sf, tol = 1e-9) {
  slo <- rep(-Inf, sf$n_sym)
  shi <- rep(Inf, sf$n_sym)
  slo[sf$var_sym] <- sf$lower
  shi[sf$var_sym] <- sf$upper
  # parameters were substituted; constants only remain
  bad <- character()
  ob <- sf$objective_lowered
  if (!is.null(ob)) {
    cv <- curvature(ob, slo, shi, tol)$curv
    want <- if (sf$sense == "min") c(CURV_CONST, CURV_AFFINE, CURV_CONVEX) else
      c(CURV_CONST, CURV_AFFINE, CURV_CONCAVE)
    if (!all(cv %in% want)) bad <- c(bad, "the objective")
  }
  off <- 0L
  for (e in sf$cons_lowered) {
    if (is.null(e)) next
    k <- e$len
    rows <- off + seq_len(k)
    cv <- curvature(e, slo, shi, tol)$curv
    up <- is.finite(sf$row_upper[rows])
    lo <- is.finite(sf$row_lower[rows])
    ok <- (!up | cv %in% c(CURV_CONST, CURV_AFFINE, CURV_CONVEX)) &
      (!lo | cv %in% c(CURV_CONST, CURV_AFFINE, CURV_CONCAVE))
    if (!all(ok)) bad <- c(bad, sf$row_names[rows][!ok])
    off <- off + k
  }
  if (!length(bad)) {
    return(list(convex = TRUE, evidence = paste(
      "structural: composition rules (disciplined convex programming) on the",
      "box of variable bounds prove the objective and every bounded side of",
      "the constraints convex")))
  }
  shown <- utils::head(bad, 5L)
  list(convex = FALSE, evidence = sprintf(paste(
    "unknown: composition rules cannot prove %s convex; solutions from local",
    "solvers are local optima"),
    paste0(paste(shown, collapse = ", "),
           if (length(bad) > 5L) sprintf(" and %d more", length(bad) - 5L))))
}
