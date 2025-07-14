# calculate W,c and the volume of B
#'
#' Determine the function W, radius c and volume of B from Metodiev et al.(2025)
#'
#' @param params      sample from the posterior (as a matrix)
#' @param limit       the limit placed on the lps
#' @param lps         values of the unnormalized log posterior density
#' @param logpost     function used to compute the lps
#' @param sims        n_simul x G x (u+1) array of parameters sampled from
#'                    the posterior, where
#'                      n_simul is the number of simulations from the posterior,
#'                      G       is the number of components,
#'                      u       is the number of mixture component parameters
#'                              (parameter u+1 is the mixture weight)
#' @param iters       half the number of simulations from the posterior
#' @param ellipse     the ellipsoid E from Metodiev et al. (2025)
#' @param type        THAMES variant ("simple", "permutations", or "standard")
#' @param lps_unif    values of the unnormalized log posterior density,
#'                    evaluated on a uniform sample on the posterior ellipsoid
#' @param max_iters   maximum number of shrinkage iterations
#' @importFrom stats          cov
#' @importFrom igraph         graph_from_adjacency_matrix shortest.paths
#' @return a named list with the following elements:
#'          perms:      number of permutations that were evaluated
#'          graph:      the overlap graph for G
#'          c_opt:      the radius of the ellipsoid in Metodiev et al. (2025)
#'          scaling:    list of fit determined via QDA (means and covariances)
#'          log_cor:    used to compute the volume of B in Metodiev et al.(2025)
#'          param_test: Monte Carlo sample used to compute the volume of B
#'          co:         the criterion of overlap for G
#'          ellipse:    the ellipsoid, potentially shifted to the sample mode
#'
#' @references  Martin Metodiev, Nicholas J. Irons, Marie Perrot-Dockès,
#' Pierre Latouche, Adrian E. Raftery. "Easily Computed Marginal Likelihoods
#' for Multivariate Mixture Models Using the THAMES Estimator."
#' arXiv preprint arXiv:2504.21812.
#'
#' @keywords internal
compute_W_c_volB = function(params,
                            limit,
                            lps,
                            logpost,
                            sims,
                            iters,
                            ellipse,
                            type="simple",
                            lps_unif=NULL,
                            max_iters=Inf){
  G = dim(sims)[2]
  theta_hat = ellipse$theta_hat
  sigma_hat = ellipse$sigma_hat
  c_opt = ellipse$c_opt
  d_par = length(theta_hat)
  graph=NULL
  co=NULL

  n_simuls = 2*iters

  simsmat = matrix(c(sims),nrow=dim(sims)[1])

  if(dim(sims)[3]==1){
    mu_post = colMeans(simsmat[1:iters,])
    sigma_post = cov(simsmat[1:iters,])
  } else{
    mu_post = colMeans(simsmat[1:iters,-ncol(simsmat)])
    sigma_post = cov(simsmat[1:iters,-ncol(simsmat)])
  }

  inv_post_var = sparsediscrim::solve_chol(sigma_post)

  counter = 0
  log_cor = -Inf

  in_ellipse = TRUE
  c_opt_old = c_opt
  center = NULL
  infinite_thames = FALSE

  # decrease radius if it takes too long to compute
  # This first loop is only to decrease computation time.
  # It works the same as the first, but does not resample via Monte Carlo.
  while(is.infinite(log_cor) & in_ellipse){
    c_opt = c_opt_old / 2^counter
    counter = counter + 1
    params_centered = params[(iters+1):(2*iters),] -
      t(matrix(rep(theta_hat,iters),ncol=iters))
    radius = c_opt

    # check if ellipse empty; set mean to mode of second half if not
    in_E =
      rowSums((params_centered %*%
                 sparsediscrim::solve_chol(ellipse$sigma_hat)) *
                params_centered) <=radius^2
    empty_ellipse =  (sum(in_E) == 0)
    log_cor = 0

    if(empty_ellipse){
      theta_hat = params[(iters+1):(2*iters),
      ][which.max(lps[(iters+1):(2*iters)]),]
      mu_post = theta_hat
      ellipse$theta_hat = theta_hat
      c_opt = sqrt(ncol(params)+1)
      log_cor = -Inf
      counter = 0
      center = t(sims[(iters+1):(2*iters),,
                      -dim(sims)[3]][which.max(lps[(iters+1):(2*iters)]),,])
    }
    if(!empty_ellipse){
      if(type=="simple"){

        # only take the sample within the ellipse if the mean had to be reset
        if(is.null(center)){
          scaling = get_qda_scaling(G,sims[(iters+1):(2*iters),,])
        } else{
          if(sum(in_E)!=1){
            scaling = get_qda_scaling(G,sims[(iters+1):(2*iters),,][in_E,,])
          } else{
            scaling = list(meanhat=center,
                           sigmahat=lapply(1:G,
                                           function(g) (1e-10)*
                                             diag(nrow(center))))
          }
        }

        graph_and_non_I_set = calc_non_I_set(scaling, G, sims, c_opt)
        scaling$non_I_set = graph_and_non_I_set$non_I_set
        if(c_opt_old == c_opt){
          overlapgraph = overlapgraph(sims)
          graph = overlapgraph$graph
          co = overlapgraph$co
        }

        # the set has to include at least two components
        if((length(scaling$non_I_set)==(G-1))){
          graph_and_non_I_set$complexity_limit_estim = factorial(G)
          scaling$non_I_set = scaling$non_I_set[-1]
        }

        # calculate an upper bound on the complexity limit using the graphmat
        graphmat_dummy = (1-graph_and_non_I_set$graphmat)

        # give all edges in the graph a direction
        graphmat_dummy[lower.tri(graphmat_dummy,diag = TRUE)] = 0

        # compute upper bound in the same way that we compute it in
        distgraph = graph_from_adjacency_matrix(graphmat_dummy)
        igraph::E(distgraph)$weight = -1
        dis <- (-shortest.paths(distgraph, v=(igraph::V(distgraph)), mode="out"))

        complexity_limit_estim = exp(lfactorial(G) - lfactorial(max(dis)))

        if((complexity_limit_estim>50000)){
          log_cor = -Inf
          if((c_opt<0)|(counter>max_iters)){
            log_cor = 0
            infinite_thames = TRUE
          }
        }
      }
    }
  }

  # exit if thames is infinite
  if(infinite_thames){
    return(list(perms=0,
                graph=graph,
                c_opt=0,
                scaling=scaling,
                log_cor=-Inf,
                param_test = params,
                co = co))
  }

  # decrease radius if none of the MCMC sample lies in the HPD region
  counter = counter - 1
  log_cor = -Inf

  in_ellipse = TRUE
  # c_opt_old = c_opt
  # center = NULL
  while(is.infinite(log_cor) & in_ellipse){
    c_opt = c_opt_old / 2^counter
    param_test = runif_in_ellipsoid(n_simuls, inv_post_var, c_opt) +
      t(matrix(rep(mu_post,n_simuls),ncol=n_simuls))
    param_test_extended = extend_param(param_test, G)
    sims_test = array(c(param_test_extended),dim=dim(sims))
    param_test = melt_sims_simple(sims_test,G)[,-(ncol(param_test)+1)]
    param_test_extended = extend_param(param_test, G)
    if((!is.null(lps_unif))&(counter==0)){
      lps_test = lps_unif
    } else{
      lps_test = logpost(sims_test)
    }

    log_cor = log(mean(lps_test>limit))
    counter = counter + 1
    params_centered = params[(iters+1):(2*iters),] -
      t(matrix(rep(theta_hat,iters),ncol=iters))
    radius = c_opt

    # check if ellipse empty; set mean to mode of second half if not
    in_E = rowSums((params_centered %*%
                      sparsediscrim::solve_chol(ellipse$sigma_hat)) *
                     params_centered) <=radius^2

    empty_ellipse =  (sum(in_E) == 0)
    if(empty_ellipse){
      theta_hat = params[(iters+1):(2*iters),
      ][which.max(lps[(iters+1):(2*iters)]),]
      mu_post = theta_hat
      ellipse$theta_hat = theta_hat
      c_opt = sqrt(ncol(params)+1)
      log_cor = -Inf
      counter = 0
      center = t(sims[(iters+1):(2*iters),,
                      -dim(sims)[3]][which.max(lps[(iters+1):(2*iters)]),,])
    }

    if(!empty_ellipse){
      if(type=="simple"){

        # only take the sample within the ellipse if the mean had to be reset
        if(is.null(center)){
          scaling = get_qda_scaling(G,sims[(iters+1):(2*iters),,])
        } else{
          if(sum(in_E)!=1){
            scaling = get_qda_scaling(G,sims[(iters+1):(2*iters),,][in_E,,])
          } else{
            scaling = list(meanhat=center,
                           sigmahat=lapply(1:G,
                                           function(g) (1e-10)*
                                             diag(nrow(center))))
          }
        }

        graph_and_non_I_set = calc_non_I_set(scaling, G, sims_test, c_opt)
        scaling$non_I_set = graph_and_non_I_set$non_I_set
        if(c_opt_old == c_opt){
          overlapgraph = overlapgraph(sims)
          graph = overlapgraph$graph
          co = overlapgraph$co
        }

        # the set has to include at least two components
        if((length(scaling$non_I_set)==(G-1))){
          graph_and_non_I_set$complexity_limit_estim = factorial(G)
          scaling$non_I_set = scaling$non_I_set[-1]
        }
        print(c_opt)
        print(graph_and_non_I_set$complexity_limit_estim)
        if((graph_and_non_I_set$complexity_limit_estim>50000)){
          log_cor = -Inf
          if((c_opt<0)|(counter>max_iters)){
            log_cor = 0
            infinite_thames = TRUE
          }
        }
      }
    }
  }

  # exit if thames is infinite
  if(infinite_thames){
    return(list(perms=0,
                graph=graph,
                c_opt=0,
                scaling=scaling,
                log_cor=-Inf,
                param_test = params,
                co = co))
  }

  if(type=="simple"){

    theta_hat_extended = extend_param(rbind(theta_hat,theta_hat),G)[1,]

    mu_post_extended = extend_param(rbind(mu_post,mu_post),G)[1,]
    sims_theta_hat_extended = array(c(mu_post_extended),dim=c(1,dim(sims)[2:3]))
    param_test_f_transform = matrix(reorder_by_lda(scaling, G, sims_test)$W,
                                    ncol=G)

    graphmat = graph_and_non_I_set$graphmat
    delta_mat = matrix(0,nrow=nrow(graphmat),ncol=ncol(graphmat))
    for(g1 in 1:(nrow(delta_mat)-1)){
      for(g2 in (g1+1):(nrow(delta_mat))){
        delta_mat[g1,g2] = (mean(param_test_f_transform[,g1] <
                                   param_test_f_transform[,g2]) == 1)
        delta_mat[g2,g1] = (mean(param_test_f_transform[,g2] <
                                   param_test_f_transform[,g1]) == 1)
      }
    }
    delta_mat = delta_mat * (1-graph_and_non_I_set$graphmat)

    adj_matrix = delta_mat
    adj_list = matrix(1:2,nrow=1)
    for(g1 in (1:nrow(adj_matrix))[rowSums(adj_matrix)>0]){
      for(g2 in which(adj_matrix[g1, ] == 1)){
        adj_list = rbind(adj_list,c(g1,g2))
      }
    }

    if(nrow(adj_list)>1){
      adj_list = lapply(2:nrow(adj_list), function(l) adj_list[l,])
    } else{
      adj_list = list()
    }

    perms = alltopsorts_recursion(G, adj_list)
    gc()
    rm()
  } else{
    perms = NULL
    scaling = NULL
  }
  return(list(perms=perms,
              graph=graph,
              c_opt=c_opt,
              scaling=scaling,
              log_cor=log_cor,
              param_test = param_test,
              co = co,
              ellipse=ellipse))
}
