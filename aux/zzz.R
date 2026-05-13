#' @keywords internal
.onLoad <- function(libname, pkgname) {
  message(sprintf("Package '%s' is being loaded from '%s'", pkgname, libname))
  required_packages <- c("GIGrvg", "statmod", "mvtnorm", "bzinb")
  missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
  if (length(missing_packages) > 0) {
    stop("Missing required packages: ", paste(missing_packages, collapse = ", "),
         ". Please install them.")
  }
}
