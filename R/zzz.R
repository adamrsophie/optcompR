#' @keywords internal
"_PACKAGE"

#' @importFrom methods as
NULL

.onLoad <- function(libname, pkgname) {
  register_highs()
  register_builtin_backends()
  invisible()
}
