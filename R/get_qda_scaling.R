# calculate the weights used for W via QDA
#'
#' Calculates the weights of the function W described in Metodiev et al. (2025).
#'
#' @param G           number of components
#' @param sims        n_simul x G x (u+1) array of parameters sampled from
#'                    the posterior, where
#'                      n_simul is the number of simulations from the posterior,
#'                      G       is the number of components,
#'                      u       is the number of mixture component parameters
#'                              (parameter u+1 is the mixture weight)
#' @param center      optional, means of the mixture component parameters
#' @return a named list including the weights (i.e., means and covariances)
#'
#' @references  Martin Metodiev, Nicholas J. Irons, Marie Perrot-Dockès,
#' Pierre Latouche, Adrian E. Raftery. "Easily Computed Marginal Likelihoods
#' for Multivariate Mixture Models Using the THAMES Estimator."
#' arXiv preprint arXiv:2504.21812.
#'
#' @keywords internal
get_qda_scaling = function(G, sims, center=NULL){
  if(length(dim(sims))==3){
    data = sapply(1:(dim(sims)[3]-1), function(r) c(sims[,,r]))
    meanhat =
      sapply(0:(G-1),
             function(g) colMeans(data[(dim(sims)[1]*g+1):(dim(sims)[1]*(g+1)),
                                       ]))
    sigmahat =
      sapply(0:(G-1),
             function(g) cov(data[(dim(sims)[1]*g+1):(dim(sims)[1]*(g+1)),]),
             simplify = FALSE)
  } else{
    sims = array(c(sims),dim=c(dim(sims),1))
    data = matrix(c(sims[,,1]),ncol=1)
    meanhat =
      sapply(0:(G-1),
             function(g) mean(data[(dim(sims)[1]*g+1):(dim(sims)[1]*(g+1)),]))
    sigmahat =
      sapply(0:(G-1),
             function(g) var(data[(dim(sims)[1]*g+1):(dim(sims)[1]*(g+1)),]),
             simplify = FALSE)
  }

  if(!is.null(center)){
    meanhat = center
  }

  scaling = list(meanhat=meanhat,sigmahat=sigmahat)
  return(scaling)
}
