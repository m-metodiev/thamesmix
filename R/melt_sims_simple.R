#' Melt array to matrix
#'
#' This function melts the simulations,
#' presented as a n_simul x G x u array,
#' into a matrix.
#'
#' @param sims n_simul x G x (u+1) array of parameters sampled from
#'                    the posterior, where
#'                      n_simul is the number of simulations from the posterior,
#'                      G       is the number of components,
#'                      u       is the number of mixture component parameters
#'                              (parameter u+1 is the mixture weight)
#' @param G   number of components
#' @return a matrix of dimension n_simul x (G(u+1)-1)
#'
#' @keywords internal
melt_sims_simple = function(sims,G){
  params = matrix(nrow=dim(sims)[1],ncol=1)
  for(index in 1:G){
    params = cbind(params,sims[,index,])
  }
  params = params[,-1]

  params = params[,rep(seq(1,G*dim(sims)[3]-(dim(sims)[3]-1),
                           dim(sims)[3]),dim(sims)[3])+
                    rep(0:(dim(sims)[3]-1),each=G)]
  return(params)
}
