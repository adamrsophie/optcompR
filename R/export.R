#' Write a model to an MPS or LP file
#'
#' Writes the standard form of a linear, quadratic or mixed-integer linear
#' model in a format other solvers read, using the HiGHS writer. The format
#' follows the file extension: `.mps` for (free) MPS or `.lp` for the CPLEX
#' LP format. Parameters are written with their current values. Variable and
#' row names are not transferred; columns and rows keep the order of
#' [opt_standard_form()].
#'
#' @param model An [opt_model()] of class LP, MILP or QP.
#' @param file Path ending in `.mps` or `.lp`.
#' @return The file path, invisibly.
#' @examples
#' m <- opt_model("max")
#' x <- opt_variable(m, "x", dim = 2, lower = 0, upper = 3)
#' opt_constraint(m, x[1] + 2 * x[2] <= 4)
#' opt_objective(m, 3 * x[1] + x[2])
#' f <- tempfile(fileext = ".mps")
#' opt_write_model(m, f)
#' head(readLines(f))
#' @export
opt_write_model <- function(model, file) {
  check_model(model)
  check_scalar_string(file, "file")
  if (!grepl("\\.(mps|lp)$", file, ignore.case = TRUE)) {
    stop("`file` must end in .mps or .lp.", call. = FALSE)
  }
  sf <- opt_standard_form(model)
  cls <- opt_classify(sf)
  if (!cls$class %in% c("LP", "MILP", "QP")) {
    stop(sprintf("Writing %s models is not supported (LP, MILP and QP only).",
                 cls$class), call. = FALSE)
  }
  types <- ifelse(sf$type %in% c("integer", "binary"), 2L, 1L)
  has_rows <- nrow(sf$A) > 0L
  Q <- if (is.null(sf$Q)) NULL else
    methods::as(methods::as(sf$Q, "generalMatrix"), "CsparseMatrix")
  hm <- highs::highs_model(
    Q = Q, L = sf$c, lower = sf$lower, upper = sf$upper,
    A = if (has_rows) methods::as(sf$A, "CsparseMatrix") else NULL,
    lhs = if (has_rows) sf$row_lower else NULL,
    rhs = if (has_rows) sf$row_upper else NULL,
    types = types, maximum = sf$sense == "max", offset = sf$offset)
  highs::highs_write_model(hm, file)
  invisible(file)
}
