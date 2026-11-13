# Text rendering of expressions.

sym_label <- function(model, j) {
  if (is.null(model)) return(paste0("z", j))
  nm <- model$sym_name[j]
  ifelse(is.na(nm), paste0("z", j), nm)
}

term_text <- function(coef, body, first) {
  if (body == "") {
    s <- fmt_num(abs(coef))
  } else if (abs(coef) == 1) {
    s <- body
  } else {
    s <- paste0(fmt_num(abs(coef)), "*", body)
  }
  if (first) {
    if (coef < 0) paste0("-", s) else s
  } else {
    paste(if (coef < 0) "-" else "+", s)
  }
}

format_poly_element <- function(e, k, drop_const = FALSE) {
  parts <- character()
  qs <- which(e$qi == k)
  for (t in qs) {
    a <- sym_label(e$model, e$qj[t])
    b <- sym_label(e$model, e$qk[t])
    body <- if (e$qj[t] == e$qk[t]) paste0(a, "^2") else paste0(a, "*", b)
    parts <- c(parts, term_text(e$qv[t], body, length(parts) == 0L))
  }
  ls <- which(e$li == k)
  for (t in ls) {
    parts <- c(parts, term_text(e$lv[t], sym_label(e$model, e$lj[t]),
                                length(parts) == 0L))
  }
  if (!drop_const && (e$const[k] != 0 || length(parts) == 0L)) {
    parts <- c(parts, term_text(e$const[k], "", length(parts) == 0L))
  }
  if (length(parts) == 0L) return("0")
  paste(parts, collapse = " ")
}

format_nl <- function(e, k) {
  if (e$type == "poly") {
    s <- format_poly_element(e, if (e$len == 1L) 1L else k)
    return(if (grepl(" ", s)) paste0("(", s, ")") else s)
  }
  arg <- function(a) format_nl(a, if (a$len == 1L) 1L else k)
  switch(e$op,
    "+" = paste0("(", arg(e$args[[1L]]), " + ", arg(e$args[[2L]]), ")"),
    "*" = paste0(arg(e$args[[1L]]), " * ", arg(e$args[[2L]])),
    "/" = paste0(arg(e$args[[1L]]), " / ", arg(e$args[[2L]])),
    "^" = paste0(arg(e$args[[1L]]), "^",
                 fmt_num(rep_len(e$data$p, e$len)[k])),
    linmap = {
      M <- e$data$M
      row <- M[k, ]
      nz <- which(row != 0)
      a <- e$args[[1L]]
      terms <- vapply(nz, function(t) {
        s <- format_nl(a, t)
        if (row[t] == 1) s else paste0(fmt_num(row[t]), "*", s)
      }, "")
      if (length(terms) == 0L) "0" else
        paste0("(", paste(terms, collapse = " + "), ")")
    },
    paste0(e$op, "(", arg(e$args[[1L]]), ")"))
}

# One line of text per element.
format_expr <- function(e, max_elements = 10L, relation = FALSE) {
  n <- min(e$len, max_elements)
  if (e$type == "poly") {
    e <- compact_poly(e)
    vapply(seq_len(n), format_poly_element, "", e = e)
  } else {
    vapply(seq_len(n), function(k) format_nl(e, k), "")
  }
}

#' @export
format.opt_expr <- function(x, ...) format_expr(x, max_elements = x$len)

#' @export
print.opt_expr <- function(x, ...) {
  shape <- if (is.null(x$dim)) as.character(x$len) else
    paste(x$dim, collapse = " x ")
  cat(sprintf("<opt_expr %s, %s>\n", shape, expr_kind_label(x)))
  txt <- format_expr(x)
  if (x$len == 1L) {
    cat(" ", txt, "\n")
  } else {
    cat(paste0("  [", seq_along(txt), "] ", txt), sep = "\n")
    if (x$len > length(txt)) cat(sprintf("  ... %d more\n", x$len - length(txt)))
  }
  invisible(x)
}

#' @export
as.character.opt_expr <- function(x, ...) format_expr(x, max_elements = x$len)
