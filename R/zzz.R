#' @keywords internal
"_PACKAGE"

#' @importFrom methods as
#' @useDynLib optcompR, .registration = TRUE
NULL

.onLoad <- function(libname, pkgname) {
  register_highs()
  register_builtin_backends()
  register_clarabel()
  register_nlpbb()
  register_ipopt()
  invisible()
}
