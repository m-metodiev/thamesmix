# calculate I(G) via quadratic programming
#'
#' Calculates the set I(G) described in Metodiev et al. (2025)
#' as well as the limit on the complexity of the THAMES associated with it.
#'
#' @param scaling     list of fit determined via QDA (means and covariances)
#' @param G           number of components
#' @param sims        n_simul x G x (u+1) array of parameters sampled from
#'                    the posterior, where
#'                      n_simul is the number of simulations from the posterior,
#'                      G       is the number of components,
#'                      u       is the number of mixture component parameters
#'                              (parameter u+1 is the mixture weight)
#' @param c_opt       radius of the ellipsoid E from Metodiev et al. (2025)
#' @importFrom stats          var
#' @importFrom sparsediscrim  solve_chol
#' @importFrom quadprog       solve.QP
#' @importFrom igraph         graph_from_adjacency_matrix shortest.paths
#' @importFrom gor            build_cover_greedy
#' @return a named list with the following elements:
#'          graph:                  the overlap graph for G
#'          non_I_set:              the complement of the set I(G)
#'          complexity_limit_estim: an upper bound on the THAMES complexity
#'          graphmat:               adjacency matrix of the overlap graph for G
#'
#' @references  Martin Metodiev, Nicholas J. Irons, Marie Perrot-Dockès,
#' Pierre Latouche, Adrian E. Raftery. "Easily Computed Marginal Likelihoods
#' for Multivariate Mixture Models Using the THAMES Estimator."
#' arXiv preprint arXiv:2504.21812.
#'
#' @keywords internal
calc_non_I_set = function(scaling, G, sims, c_opt){
  meanhat = scaling$meanhat
  sigmahat = scaling$sigmahat
  if(dim(sims)[3]>1){
    data = sapply(1:(dim(sims)[3]-1), function(r) c(sims[,,r]))

    simsmat = matrix(c(sims),nrow=dim(sims)[1])
    mu_post = colMeans(simsmat[1:(dim(sims)[1]/2),-ncol(simsmat)])
    sigma_post_inverse =
      sparsediscrim::solve_chol(cov(simsmat[1:(dim(sims)[1]/2),-ncol(simsmat)]))
  } else{
    data = matrix(c(sims[,,1]),ncol=1)
    meanhat =
      sapply(0:(G-1),
             function(g) mean(data[(dim(sims)[1]*g+1):(dim(sims)[1]*(g+1)),]))
    sigmahat =
      sapply(0:(G-1),
             function(g) var(data[(dim(sims)[1]*g+1):(dim(sims)[1]*(g+1)),]),
             simplify = FALSE)

    simsmat = matrix(c(sims),nrow=dim(sims)[1])
    mu_post = colMeans(simsmat[1:(dim(sims)[1]/2),])
    sigma_post_inverse =
      sparsediscrim::solve_chol(cov(simsmat[1:(dim(sims)[1]/2),]))
  }

  graphmat = diag(G)

  for(g1 in 1:(G-1)){
    for(g2 in (g1+1):G){
      Amat = matrix(0,nrow=nrow(sigma_post_inverse),
                    ncol=ncol(sigma_post_inverse))
      Amat_func = function(u){
        sims_identification = array(0,dim=c(1,dim(sims)[-1]))
        sims_identification[1,g1,u] = 1
        sims_identification[1,g2,u] = -1
        component_pos = matrix(c(sims_identification),nrow=1)[1,]
        if(dim(sims)[3]==1){
          return(component_pos)
        } else{
          return(component_pos[-length(component_pos)])
        }
      }

      if(dim(sims)[3]==1){
        Amat_eqs = t(simplify2array(lapply(1:dim(sims)[3],Amat_func)))
      } else{
        Amat_eqs = t(simplify2array(lapply(1:(dim(sims)[3]-1),Amat_func)))
      }

      Amat[1:nrow(Amat_eqs),] = Amat_eqs

      divide_to_stabilize = 1
      solution_of_QP =
        try(quadprog::solve.QP(Dmat=sigma_post_inverse,
                               dvec=sigma_post_inverse%*%mu_post,
                               Amat = t(Amat),
                               bvec = rep(0,ncol(Amat)),
                               meq = ncol(Amat)))
      # solving some numerical issues
      while(is.character(solution_of_QP)){
        divide_to_stabilize = divide_to_stabilize * 2
        solution_of_QP =
          try(quadprog::solve.QP(Dmat=sigma_post_inverse/divide_to_stabilize,
                                 dvec=(sigma_post_inverse%*%mu_post)/
                                   divide_to_stabilize,
                                 Amat = t(Amat),
                                 bvec = rep(0,ncol(Amat)),
                                 meq = ncol(Amat)))
      }

      maxnorm_ellipse = 2*solution_of_QP$value * divide_to_stabilize +
        t(mu_post)%*%sigma_post_inverse%*%t(t(mu_post))
      graphmat[g1,g2] = 0 + (maxnorm_ellipse <= (c_opt^2))
      graphmat[g2,g1] = 0 + (maxnorm_ellipse <= (c_opt^2))
    }
  }

  diag(graphmat)=0
  graph = graph_from_adjacency_matrix(graphmat,mode="undirected")
  non_I_set = gor::build_cover_greedy(graph)$set
  igraph::V(graph)$color = rep("white",G)
  igraph::V(graph)$label = rep(" ",G)

  Wmat =
    simplify2array(lapply(1:G,
                          function(g) param_i_qda_linearized(g, sims,
                                                             meanhat,
                                                             sigmahat,
                                                             non_I_set)))
  W=c(Wmat)

  ranges = sapply(0:(G-1),
                  function(i) range(W[(i*dim(sims)[1]+1):((i+1)*dim(sims)[1])]))

  delta_mat = matrix(0,nrow=nrow(graphmat),ncol=ncol(graphmat))
  for(g1 in 1:(nrow(delta_mat)-1)){
    for(g2 in (g1+1):(nrow(delta_mat))){
      delta_mat[g1,g2] = (ranges[2,g1] < ranges[1,g2]) + 0
    }
  }

  distgraph = graph_from_adjacency_matrix(delta_mat)
  igraph::E(distgraph)$weight = -1
  dis <- (-shortest.paths(distgraph, v=(igraph::V(distgraph)), mode="out"))

  complexity_limit_estim = exp(lfactorial(G) - lfactorial(max(dis)))
  return(list(graph=graph,
              non_I_set=non_I_set,
              complexity_limit_estim=complexity_limit_estim,
              graphmat=graphmat))
}
