# Optional interface to symcompR symbolic expressions.
#
# symcompR expressions are lists of class "sym" with fields `type`, `args`,
# `name` and `value` (Add, Mul, Pow, Symbol, Integer, Rational, Float,
# Function nodes). Conversion is exact for these node types; numbers are
# converted to double precision.

need_symcompR <- function() {
  if (!requireNamespace("symcompR", quietly = TRUE)) {
    stop("This function needs the symcompR package (remotes::install_github('adamrsophie/symcompR')).",
         call. = FALSE)
  }
}

sym_number <- function(node) {
  v <- node$value
  if (node$type %in% c("Integer", "Rational")) return(as.numeric(v))
  if (node$type == "Float") return(as.numeric(v))
  NULL
}

#' Convert a symcompR expression to an optimization expression
#'
#' Translates a symbolic expression from the symcompR package into an
#' `opt_expr` of a model. Symbols are matched by name to the model's
#' variables and parameters (scalar elements such as `x[2]` included), or to
#' the expressions in `map`. Supported nodes are sums, products, powers,
#' integers, rationals, floats, numeric constants such as pi, and the
#' functions `exp`, `log`, `sqrt`, `sin`, `cos`, `tan`, `tanh`, `sinh`,
#' `cosh`, `asin`, `acos`, `atan` and `Abs`. Any other node is an error, so
#' a converted expression always means the same function as the symbolic
#' one. Rational numbers become double-precision values.
#'
#' @param model An [opt_model()].
#' @param expr A symcompR expression (class `sym`) or number.
#' @param map Optional named list from symbol names to expressions or
#'   numbers; it takes precedence over model names.
#' @return An `opt_expr`.
#' @seealso [opt_to_sym()]
#' @examplesIf requireNamespace("symcompR", quietly = TRUE)
#' library(symcompR)
#' m <- opt_model()
#' x <- opt_variable(m, "x")
#' y <- opt_variable(m, "y")
#' s <- (Symbol("x") - 1)^2 + Symbol("x") * Symbol("y") + exp(Symbol("y"))
#' e <- opt_from_sym(m, s)
#' opt_objective(m, e)
#' opt_classify(m)$class
#' @export
opt_from_sym <- function(model, expr, map = list()) {
  check_model(model)
  need_symcompR()
  if (is.numeric(expr)) return(as_expr(expr))
  if (!inherits(expr, "sym")) stop("`expr` must be a symcompR expression.",
                                   call. = FALSE)
  out <- sym_convert(expr, model, map)
  if (is_expr(out) && is.null(out$model)) out$model <- NULL
  out
}

sym_convert <- function(node, model, map) {
  num <- sym_number(node)
  if (!is.null(num)) return(as_expr(num))
  ty <- node$type
  args <- node$args
  if (ty == "Symbol") {
    nm <- node$name
    if (!is.null(map[[nm]])) return(as_expr(map[[nm]]))
    j <- match(nm, model$sym_name)
    if (is.na(j)) {
      bid <- index_get(model$block_index, nm)
      if (is.na(bid) || length(model$blocks[[bid]]$sym) != 1L) {
        stop(sprintf("Symbol '%s' is not a scalar variable or parameter of the model; supply it in `map`.",
                     nm), call. = FALSE)
      }
      j <- model$blocks[[bid]]$sym
    }
    return(symbol_expr(model, j, NULL))
  }
  conv <- lapply(args, sym_convert, model = model, map = map)
  if (ty == "Add") return(Reduce(expr_add, conv))
  if (ty == "Mul") return(Reduce(expr_mul, conv))
  if (ty == "Pow") return(expr_pow(conv[[1L]], conv[[2L]]))
  if (ty == "Function") {
    fn <- node$name
    fn <- if (identical(fn, "Abs")) "abs" else fn
    if (fn %in% unary_ops && length(conv) == 1L) return(expr_unary(fn, conv[[1L]]))
    stop(sprintf("The symbolic function '%s' has no optcompR counterpart.", fn),
         call. = FALSE)
  }
  if (length(args) == 0L) {
    v <- tryCatch(as.numeric(symcompR::N(node)$value), error = function(e) NULL)
    if (length(v) == 1L && is.finite(v)) return(as_expr(v))
  }
  stop(sprintf("Cannot convert symbolic node of type '%s'.", ty), call. = FALSE)
}

#' Convert an optimization expression to a symcompR expression
#'
#' Builds the symbolic form of a scalar expression, with one symcompR
#' `Symbol` per variable or parameter (named as in the model). Coefficients
#' are double-precision floats unless they are integers. The symbolic form
#' can be simplified or differentiated with symcompR and converted back with
#' [opt_from_sym()].
#'
#' @param x A scalar `opt_expr`.
#' @return A symcompR expression.
#' @examplesIf requireNamespace("symcompR", quietly = TRUE)
#' m <- opt_model()
#' x <- opt_variable(m, "x")
#' s <- opt_to_sym(3 * x^2 + exp(x))
#' diff(s, symcompR::Symbol("x"))
#' @export
opt_to_sym <- function(x) {
  need_symcompR()
  x <- as_expr(x)
  if (x$len != 1L) stop("opt_to_sym() needs a scalar expression.",
                        call. = FALSE)
  to_sym_node(x, 1L)
}

sym_num <- function(v) {
  if (v == round(v) && abs(v) < 2^31) symcompR::Integer(as.integer(v)) else
    symcompR::Float(v)
}

to_sym_node <- function(e, k) {
  S <- function(j) symcompR::Symbol(e$model$sym_name[j])
  if (e$type == "poly") {
    e <- compact_poly(e)
    out <- sym_num(e$const[k])
    for (t in which(e$li == k)) out <- out + sym_num(e$lv[t]) * S(e$lj[t])
    for (t in which(e$qi == k)) {
      out <- out + sym_num(e$qv[t]) * S(e$qj[t]) * S(e$qk[t])
    }
    return(out)
  }
  arg <- function(a) to_sym_node(a, if (a$len == 1L) 1L else k)
  switch(e$op,
    "+" = arg(e$args[[1L]]) + arg(e$args[[2L]]),
    "*" = arg(e$args[[1L]]) * arg(e$args[[2L]]),
    "/" = arg(e$args[[1L]]) / arg(e$args[[2L]]),
    "^" = arg(e$args[[1L]])^sym_num(rep_len(e$data$p, e$len)[k]),
    linmap = {
      row <- e$data$M[k, ]
      out <- sym_num(0)
      for (t in which(row != 0)) {
        out <- out + sym_num(row[t]) * to_sym_node(e$args[[1L]], t)
      }
      out
    },
    fn = stop(sprintf("External function '%s' has no symbolic form.",
                      e$data$name), call. = FALSE),
    # base math functions dispatch to symcompR's Math method
    match.fun(e$op)(arg(e$args[[1L]])))
}
