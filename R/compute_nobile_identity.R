#' Nobile's identity for the marginal likelihood
#'
#' This function uses the identity from Nobile (2004, 2007) to compute an
#' estimate of the marginal likelihood for a mixture model with G components
#' given an estimate of the marginal likelihood for a mixture model with G-1
#' components and an estimate of the proportion of empty components.
#'
#' @param logZhatGminus1  estimate of the marginal likelihood for G-1
#' @param p0hat_value    estimate of the proportion of empty components
#' @param G               number of components
#' @param dirichlet_vec   hyperparameter-vector of the dirichlet prior
#' @param n               size of the data
#'
#' @return estimate of the marginal likelihood for G
#'
#'@examples
#' # computes log marginal likelihood of the Swiss banknote dataset
#' # for G=4, given the settings in Metodiev et al. (2025)
#' compute_nobile_identity(logZhatGminus1 = -909.49,
#' p0hat_value = 1/4,
#' dirichlet_vec = rep(1,4),
#' n=200)
#'
#' @references Nobile, A. (2004).  On the posterior distribution of the number
#' of components in a finite mixture. The Annals of Statistics 32(5), 2044–2073.
#'
#' Nobile, A. (2007). Bayesian finite mixtures: a note on prior specification
#' and posterior computation.arXiv preprint arXiv:0711.0458.
#'
#' Martin Metodiev, Nicholas J. Irons, Marie Perrot-Dockès,
#' Pierre Latouche, Adrian E. Raftery. "Easily Computed Marginal Likelihoods
#' for Multivariate Mixture Models Using the THAMES Estimator."
#' arXiv preprint arXiv:2504.21812.
#'
#' @export
compute_nobile_identity =
  function(logZhatGminus1,p0hat_value,G,dirichlet_vec,n){
    return(lgamma(sum(dirichlet_vec))+
           lgamma(n+sum(dirichlet_vec[length(dirichlet_vec)]))-
           lgamma(sum(dirichlet_vec[-length(dirichlet_vec)]))-
           lgamma(n+sum(dirichlet_vec))+
           logZhatGminus1-
           log(p0hat_value))
}
