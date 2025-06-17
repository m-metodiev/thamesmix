# calculate W via QDA for all component parameters
#'
#' Calculates the function W described in Metodiev et al. (2025).
#'
#' @param scaling     list of fit determined via QDA (means and covariances)
#' @param G           number of components
#' @param sims        n_simul x G x (u+1) array of parameters sampled from
#'                    the posterior, where
#'                      n_simul is the number of simulations from the posterior,
#'                      G       is the number of components,
#'                      u       is the number of mixture component parameters
#'                              (parameter u+1 is the mixture weight)
#' @return a named list with the values of W as its element
#'
#' @references  Martin Metodiev, Nicholas J. Irons, Marie Perrot-Dockès,
#' Pierre Latouche, Adrian E. Raftery. "Easily Computed Marginal Likelihoods
#' for Multivariate Mixture Models Using the THAMES Estimator."
#' arXiv preprint arXiv:2504.21812.
#'
#' @keywords internal
reorder_by_lda = function(scaling, G, sims){
  meanhat = scaling$meanhat
  sigmahat = scaling$sigmahat
  non_I_set = scaling$non_I_set
  if((length(dim(sims))==3) & (dim(sims)[3]>1)){
    data = sapply(1:(dim(sims)[3]-1), function(r) c(sims[,,r]))
  } else{
    data = matrix(c(array(sims,dim=c(dim(sims)[1:2],1))[,,1]),ncol=1)
  }

  Wmat = simplify2array(lapply(1:G,
                               function(g) param_i_qda_linearized(g, sims,
                                                                  meanhat,
                                                                  sigmahat,
                                                                  non_I_set)))
  W=c(Wmat)
  return(list(W=W))
}
