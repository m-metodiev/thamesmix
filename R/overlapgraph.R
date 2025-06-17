#' Estimator of the overlap graph
#'
#' This function computes the overlap graph for mixture models.
#'
#' @param sims        n_simul x G x (u+1) array of parameters sampled from
#'                    the posterior, where
#'                      n_simul is the number of simulations from the posterior,
#'                      G       is the number of components,
#'                      u       is the number of mixture component parameters
#'                              (parameter u+1 is the mixture weight)
#'
#' @return Returns a named list with the following elements:
#'
#'          graph,              the overlap graph
#'
#'          co,                 the criterion of overlap
#'
#' @examples
#' # toy sample from the posterior
#' mus = rbind(c(17.67849, 21.46734),
#'             c(17.67849, 21.46734),
#'             c(16.98067, 21.11391),
#'             c(20.58628, 21.22104),
#'             c(17.38332, 21.37224),
#'             c(16.43644, 21.19085),
#'             c(19.49676, 21.28964),
#'             c(17.82287, 21.22475),
#'             c(18.06050, 21.36945),
#'             c(18.70759, 21.60244),
#'             c(15.93795, 21.04681),
#'             c(16.23184, 20.96049))
#' sigmasqus = rbind(c(46.75089, 3.660171),
#'                   c(58.44208, 3.026577),
#'                   c(63.19334, 4.090872),
#'                   c(87.02758, 2.856063),
#'                   c(82.34268, 3.760550),
#'                   c(50.92386, 2.380784),
#'                   c(49.51412, 3.605798),
#'                   c(38.67681, 3.362407),
#'                   c(49.59170, 3.130254),
#'                   c(63.41569, 2.475669),
#'                   c(65.95225, 3.927501),
#'                   c(47.22989, 5.465702))
#' taus = rbind(c(0.2653882, 0.7346118),
#'              c(0.2560075, 0.7439925),
#'              c(0.2371868, 0.7628132),
#'              c(0.2998265, 0.7001735),
#'              c(0.3518301, 0.6481699),
#'              c(0.2840316, 0.7159684),
#'              c(0.2060193, 0.7939807),
#'              c(0.2859257, 0.7140743),
#'              c(0.2420695, 0.7579305),
#'              c(0.2466622, 0.7533378),
#'              c(0.2726186, 0.7273814),
#'              c(0.2738916, 0.7261084))
#' sims = array(dim=c(12,2,3))
#' sims[,,1] = mus
#' sims[,,2] = sigmasqus
#' sims[,,3] = taus
#'
#' overlapgraph(sims)$co
#'
#' @references  Martin Metodiev, Nicholas J. Irons, Marie Perrot-Dockès,
#' Pierre Latouche, Adrian E. Raftery. "Easily Computed Marginal Likelihoods
#' for Multivariate Mixture Models Using the THAMES Estimator."
#' arXiv preprint arXiv:2504.21812.
#'
#' @export
overlapgraph = function(sims){
  iters = floor(dim(sims)[1]/2)
  G = dim(sims)[2]
  d_par = dim(sims)[2]*dim(sims)[3] - 1
  scaling = get_qda_scaling(G,sims[(iters+1):(2*iters),,],center=NULL)
  graph_and_non_I_set = calc_non_I_set(scaling, G, sims, c_opt=sqrt(d_par+1))

  scaling$non_I_set = graph_and_non_I_set$non_I_set
  len_I_set = G - length(scaling$non_I_set)

  return(list(graph=graph_and_non_I_set$graph,
              co = len_I_set - length(scaling$non_I_set)))
}
