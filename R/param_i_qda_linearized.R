# calculate W via QDA for one component parameter
#'
#' Calculates the function W described in Metodiev et al. (2025).
#'
#' @param g           number indicating the component g
#' @param sims        n_simul x G x (u+1) array of parameters sampled from
#'                    the posterior, where
#'                      n_simul is the number of simulations from the posterior,
#'                      G       is the number of components,
#'                      u       is the number of mixture component parameters
#'                              (parameter u+1 is the mixture weight)
#' @param meanhat     list of means of component parameters
#' @param sigmahat    list of covariance matrices of component parameters
#' @param non_I_set   complement of an independent set of a DAG
#' @importFrom stats    dnorm
#' @importFrom mvtnorm  dmvnorm
#' @return a vector: W evaluated at the posterior sample of component g
#'
#' @references  Martin Metodiev, Nicholas J. Irons, Marie Perrot-Dockès,
#' Pierre Latouche, Adrian E. Raftery. "Easily Computed Marginal Likelihoods
#' for Multivariate Mixture Models Using the THAMES Estimator."
#' arXiv preprint arXiv:2504.21812.
#'
#' @keywords internal
param_i_qda_linearized = function(g, sims, meanhat, sigmahat, non_I_set){
  G = length(sigmahat)
  if((length(dim(sims))==3) & (dim(sims)[3]>1)){
    testmat = sims[,g,-dim(sims)[3]]
    postprob_mat = sapply(1:G,
                          function(g_2) mvtnorm::dmvnorm(testmat,
                                                         mean=meanhat[,g_2],
                                                         sigma=sigmahat[[g_2]],
                                                         log=TRUE))
  } else{
    testvec = array(sims,dim=c(dim(sims)[1:2],1))[,g,]
    postprob_mat = sapply(1:G,
                          function(g_2) dnorm(testvec,
                                              mean=meanhat[g_2],
                                              sd=sqrt(sigmahat[[g_2]]),
                                              log=TRUE))
  }

  if(is.matrix(postprob_mat)){
    postprob_mat[,non_I_set] = -Inf
    # normalize to deal with potential numeric issues
    maxlogrows = do.call(pmax, c(as.data.frame(postprob_mat)))

    postprob_mat = postprob_mat - maxlogrows
    postprob_mat_normalized = exp(postprob_mat) / rowSums(exp(postprob_mat))

    # taking the argmax (linearized)
    helper_mat = t(matrix(rep(1:G,nrow(postprob_mat)),nrow=G)) *
      (postprob_mat==0)
    g_hat = do.call(pmax, c(as.data.frame(helper_mat)))

    maxlogrowsnormalized = do.call(pmax,
                                   c(as.data.frame(postprob_mat_normalized)))
    Wmat_row_g = (g_hat + (1-maxlogrowsnormalized))

  } else{
    postprob_mat[non_I_set] = -Inf
    maxlogrows = max(postprob_mat)
    postprob_mat = postprob_mat - maxlogrows
    postprob_mat_normalized = exp(postprob_mat) / sum(exp(postprob_mat))

    g_hat = which.max(postprob_mat)

    maxlogrowsnormalized = max(postprob_mat_normalized)
    Wmat_row_g = (g_hat + (1-maxlogrowsnormalized))
  }
  return(Wmat_row_g)

}
