#' Default half-Cauchy scale for the TPB global variance parameter
#'
#' This uses the Piironen-Vehtari global-shrinkage scaling idea with a fixed
#' sparsity guess p0 = 0.01D, so p0 / (D - p0) = 0.01 / 0.99. The response is
#' not standardized internally; when y is supplied its sample sd enters the
#' scale. The final value is capped at 1.
#'
#' @param X Optional design matrix.
#' @param n Optional sample size, used when X is not supplied.
#' @param y Optional response vector. If supplied, sd(y) is used as the scale
#'   multiplier without standardizing y.
#' @return A positive numeric scale for phi.
#' @export
tpb_default_scale_phi <- function(X = NULL, y = NULL, n = NULL) {
  if (!is.null(X)) {
    n <- nrow(X)
  }
  if (is.null(n) || !is.finite(n) || n <= 0) {
    return(1e-2)
  }
  sigma_reference <- if (!is.null(y)) {
    stats::sd(as.numeric(y), na.rm = TRUE)
  } else {
    1
  }
  if (!is.finite(sigma_reference) || sigma_reference <= 0) {
    sigma_reference <- 1
  }

  odds <- 0.01 / 0.99
  pmin(odds * sigma_reference / sqrt(n), 1)
}
