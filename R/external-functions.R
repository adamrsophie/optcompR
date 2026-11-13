# External scalar functions of expressions: f(x) given as R code.

fd_gradient <- function(f, v, h) {
  vapply(seq_along(v), function(i) {
    step <- h * max(1, abs(v[i]))
    e <- replace(numeric(length(v)), i, step)
    (f(v + e) - f(v - e)) / (2 * step)
  }, 0)
}

#' Wrap an R function as a nonlinear expression
#'
#' Creates a scalar expression \eqn{f(x)} from an R function of a numeric
#' vector, for objectives or constraints that are not written with optcompR
#' operators: simulation-based losses, functions written for jaxR, or
#' functions of a model from modelcompR. Solvers need a gradient, and its
#' source is explicit:
#'
#' * a function returning the gradient vector (exact, as supplied);
#' * `"jaxR"`: reverse-mode automatic differentiation with
#'   `jaxR::grad()`. `fn` must be traceable by jaxR and receives a jaxR
#'   array, whose indexing is 0-based as in JAX (`v[0]` is the first
#'   element), so the same `fn` is evaluated by jaxR for values too. It
#'   runs with 64-bit types enabled (`jax_enable_x64`). jaxR has no
#'   prebuilt XLA runtime for Windows;
#' * `"finite_difference"`: central differences with relative step `step`.
#'   This is an approximation with error of order `step^2` (plus rounding
#'   of order `eps / step`); results record it, and it is never chosen
#'   implicitly.
#'
#' Hessians are available only when `hessian` is a function.
#'
#' @param fn Function of a numeric vector returning a number.
#' @param x Expression the function is applied to (for example a variable
#'   vector).
#' @param gradient A gradient function, `"jaxR"` or `"finite_difference"`.
#' @param hessian Optional function returning the Hessian matrix.
#' @param name Label used in printing.
#' @param step Relative step for finite differences.
#' @return A scalar nonlinear `opt_expr`.
#' @examples
#' m <- opt_model()
#' x <- opt_variable(m, "x", dim = 2, start = c(-1, 1))
#' rosen <- function(v) (1 - v[1])^2 + 100 * (v[2] - v[1]^2)^2
#' rosen_grad <- function(v) c(-2 * (1 - v[1]) - 400 * v[1] * (v[2] - v[1]^2),
#'                             200 * (v[2] - v[1]^2))
#' opt_objective(m, opt_function(rosen, x, gradient = rosen_grad,
#'                               name = "rosenbrock"))
#' if (requireNamespace("nloptr", quietly = TRUE)) opt_value(opt_solve(m), x)
#' @export
opt_function <- function(fn, x, gradient, hessian = NULL, name = "f",
                         step = 1e-6) {
  if (!is.function(fn)) stop("`fn` must be a function.", call. = FALSE)
  x <- as_expr(x)
  check_scalar_string(name, "name")
  source <- if (is.function(gradient)) "supplied" else
    match.arg(gradient, c("jaxR", "finite_difference"))
  g <- switch(source,
    supplied = gradient,
    jaxR = {
      if (!requireNamespace("jaxR", quietly = TRUE)) {
        stop("gradient = \"jaxR\" needs the jaxR package.", call. = FALSE)
      }
      jg <- jaxR::grad(fn)
      function(v) jax64(as.numeric(jg(jaxR::as_jax(v))))
    },
    finite_difference = function(v) fd_gradient(fn, v, step))
  if (!is.null(hessian) && !is.function(hessian)) {
    stop("`hessian` must be a function or NULL.", call. = FALSE)
  }
  raw <- if (source == "jaxR") function(v) jax64(as.numeric(fn(jaxR::as_jax(v))))
    else fn
  value <- function(v) {
    out <- raw(v)
    if (length(out) != 1L || !is.finite(out)) {
      stop(sprintf("External function '%s' must return one finite number.",
                   name), call. = FALSE)
    }
    out
  }
  new_nl(x$model, 1L, NULL, "fn", list(x),
         list(f = value, grad = g, hess = hessian, name = name,
              source = source))
}

fn_node_eval <- function(e, v) e$data$f(v)

# jaxR arrays are 32-bit by default, as in JAX; optimization needs doubles,
# so jaxR code runs with 64-bit types enabled and the setting is restored.
jax64 <- function(code) jaxR::with_jax_config("jax_enable_x64", TRUE, code)

fn_node_grad <- function(e, v) {
  g <- as.numeric(e$data$grad(v))
  if (length(g) != length(v)) {
    stop(sprintf("Gradient of '%s' has length %d, expected %d.", e$data$name,
                 length(g), length(v)), call. = FALSE)
  }
  g
}

# Nodes using finite differences, for labeling results.
expr_fd_nodes <- function(e) {
  if (e$type == "poly") return(character())
  own <- if (identical(e$op, "fn") && e$data$source == "finite_difference")
    e$data$name else character()
  unique(c(own, unlist(lapply(e$args, expr_fd_nodes))))
}
