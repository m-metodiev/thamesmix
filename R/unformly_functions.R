#' Uniform sampling in ellipsoid
#' @description Uniform sampling on an ellipsoid or in an ellipsoid.
#'   The sampling \emph{in} an ellipsoid is available in arbitrary
#'   dimension. The sampling \emph{on} an ellipsoid is available only in
#'   dimension 2 or 3.
#'
#' @rdname runif_ellipsoid
#' @keywords internal
runif_in_ellipsoid <- function(n, A, r){
  U <- chol(A)
  x <- runif_in_sphere(n, d=ncol(A), r=r)
  t(backsolve(U, t(x)))
}

#' Uniform sampling on/in sphere
#' @description Uniform sampling on a sphere or in a sphere, in arbitrary
#' dimension.
#' @name runif_sphere
#' @rdname runif_sphere
#'
#' @param n number of simulations
#' @param d dimension of the space
#' @param r radius of the sphere
#'
#' @return The simulations in a \code{n} times \code{d} matrix.
#' @importFrom stats runif rnorm
#'
NULL

#' @rdname runif_sphere
#' @export
runif_on_sphere <- function(n, d, r=1){
  sims <- matrix(rnorm(n*d), nrow=n, ncol=d)
  r * sims / sqrt(apply(sims, 1L, crossprod))
}

#' @rdname runif_sphere
#' @keywords internal
runif_in_sphere <- function(n, d, r=1){
  r * runif_on_sphere(n, d, runif(n)^(1/d))
}
