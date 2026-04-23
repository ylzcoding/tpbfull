#' @keywords internal
.onLoad <- function(libname, pkgname) {
  message(sprintf("Package '%s' is being loaded from '%s'", pkgname, libname))
  required_packages <- c("GIGrvg", "statmod", "mvtnorm", "coda", "gsl")
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
  ]
  if (length(missing_packages) > 0) {
    stop("Missing required packages: ", paste(missing_packages, collapse = ", "),
         ". Please install them.")
  }
}
