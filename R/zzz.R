#' @keywords internal
.onLoad <- function(libname, pkgname) {
  required_packages <- c("statmod", "mvtnorm", "coda", "gsl", "GIGrvg", "natural")
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
  ]
  if (length(missing_packages) > 0) {
    stop("Missing required packages: ", paste(missing_packages, collapse = ", "),
         ". Please install them.")
  }
}
