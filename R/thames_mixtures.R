#' THAMES estimator of the reciprocal log marginal likelihood for mixture models
#'
#' This function computes the THAMES estimate of the reciprocal
#'     log marginal likelihood for mixture models using posterior samples and
#'     unnormalized log posterior values.
#'
#' @param logpost     function logpost(sims,G) to compute lps with input "sims"
#' @param sims        n_simul x G x (u+1) array of parameters sampled from
#'                    the posterior, where
#'                      n_simul is the number of simulations from the posterior,
#'                      G       is the number of components,
#'                      u       is the number of mixture component parameters
#'                              (parameter u+1 is the mixture weight)
#' @param n_samples   integer, number of posterior samples
#' @param c_opt       radius of the ellipsoid used to compute the THAMES
#' @param type        THAMES variant ("simple", "permutations", or "standard")
#' @param seed        a seed
#' @param lps         values of the unnormalized log posterior density
#' @param lps_unif    values of the unnormalized log posterior density,
#'                    evaluated on a uniform sample on the posterior ellipsoid
#' @param max_iters   maximum number of shrinkage iterations
#' @importFrom Rfast        rowOrder rep_row
#' @importFrom stats        median
#' @importFrom combinat permn
#' @importFrom withr with_seed
#'
#' @examples
#'
#' y = c(9.172, 9.350, 9.483, 9.558, 9.775, 10.227, 10.406, 16.084, 16.170,
#'   18.419, 18.552, 18.600, 18.927, 19.052, 19.070, 19.330, 19.343, 19.349,
#'   19.440, 19.473, 19.529, 19.541, 19.547, 19.663, 19.846, 19.856, 19.863,
#'   19.914, 19.918, 19.973, 19.989, 20.166, 20.175, 20.179, 20.196, 20.215,
#'   20.221, 20.415, 20.629, 20.795, 20.821, 20.846, 20.875, 20.986, 21.137,
#'   21.492, 21.701, 21.814, 21.921, 21.960, 22.185, 22.209, 22.242, 22.249,
#'   22.314, 22.374, 22.495, 22.746, 22.747, 22.888, 22.914, 23.206, 23.241,
#'   23.263, 23.484, 23.538, 23.542, 23.666, 23.706, 23.711, 24.129, 24.285,
#'   24.289, 24.366, 24.717, 24.990, 25.633, 26.690, 26.995, 32.065, 32.789,
#'   34.279)
#'
#' R <- diff(range(y))
#' m <- mean(range(y))
#'
#' # likelihood
#' loglik_gmm <- function(sims,G){
#'   mus = sims[,,1]
#'   sigma_squs = sims[,,2]
#'   pis = sims[,,3]
#'   log_single_y = Vectorize(function(x)
#'     log(rowSums(sapply(1:G,
#'       function(g) pis[,g]*dnorm(x,mus[,g],sqrt(sigma_squs[,g]))))
#'     )
#'   )
#'   res = suppressWarnings(rowSums(log_single_y(y)))
#'   return(rowSums(log_single_y(y)))
#' }
#'
#' # prior
#' logprior_gmm_marginal <- function(sims,G) {
#'   mus = sims[,,1]
#'   sigma_squs = sims[,,2]
#'   pis = sims[,,3]
#'
#'   l_mus <- rowSums(sapply(1:G, function(g) dnorm(mus[,g], mean = m, sd = R,
#'                                                   log = TRUE)))
#'   l_pis <- LaplacesDemon::ddirichlet(1:G/G, rep(1,G),log=TRUE)
#'   l_sigma_squs <- lgamma(2*G+0.2) - lgamma(0.2) +
#'     0.2*log(10/R^2) - (2*G+0.2) * log(rowSums(sigma_squs^(-1))+10/R^2) -
#'     3*rowSums(log(sigma_squs))
#'   return(l_mus + l_pis + l_sigma_squs)
#' }
#' # unnormalized log-posterior density
#' logpost = function(sims){
#'   G = dim(sims)[2]
#'   mus = sims[,1:G,1]
#'   # apply exp transform
#'   sims[,1:G,2] = sims[,1:G,2]
#'   sigma_squs = sims[,1:G,2]
#'   pis = sims[,1:G,3]
#'
#'   # set to 0 outside of support
#'   if(G>2){
#'     mask = (((pis > 0) & (rowSums(pis[,1:(G-1)])<=1)) & (sigma_squs>0))
#'   }else{
#'     mask = (((pis > 0) & (pis[,1]<=1)) & (sigma_squs>0))
#'   }
#'   l_total = suppressWarnings(loglik_gmm(sims,G)+
#'     logprior_gmm_marginal(sims,G))
#'   l_total[exp(rowSums(log(mask)))==0] = -Inf
#'   return(l_total)
#' }
#'
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
#' # estimate of the log marginal likelihood
#' -thames_mixtures(logpost,sims)$log_zhat_inv
#'
#'
#'
#' @return Returns a named list with the following elements:
#'
#'          theta_hat,          posterior mean
#'
#'          sigma_hat,          posterior covariance matrix
#'
#'          log_det_sigma_hat,  log-determinant of sigma_hat
#'
#'          logvolA,            log-volume of the ellipsoid
#'
#'          log_zhat_inv,       log-reciprocal-marginal likelihood
#'
#'          log_zhat_inv_L,     lower bound
#'
#'          log_zhat_inv_U,     upper bound
#'
#'          alpha,              HPD-region correction
#'
#'          len_perms,          number of permutations evaluated
#'
#'          log_cor,            log-correction of the volume of the ellipsoid
#'
#'          etas,               Monte-Carlo sample on the ellipsoid
#'
#'          graph,              the overlap graph for G
#'
#'          se,                 standard_error
#'
#'          phi,                ar(1) model parameter
#'
#'          c_opt,              radius of the ellipsoid
#'
#'          d_par,              dimension of the parameter
#'
#'          G,                  number of mixture components
#'
#'          scaling,            list of fit of QDA (means, covariances)
#'
#'          co,                 the criterion of overlap
#'
#' @references  Martin Metodiev, Nicholas J. Irons, Marie Perrot-Dockès,
#' Pierre Latouche, Adrian E. Raftery. "Easily Computed Marginal Likelihoods
#' for Multivariate Mixture Models Using the THAMES Estimator."
#' arXiv preprint arXiv:2504.21812.
#'
#' @export
thames_mixtures <- function(logpost,
                            sims,
                            n_samples = NULL,
                            c_opt = NULL,
                            type="simple",
                            seed=NULL,
                            lps=NULL,
                            lps_unif=NULL,
                            max_iters=Inf
){
  if(!is.null(n_samples)){
    sims = sims[n_samples,,]
  }

  if(is.null(lps)){
    lps = logpost(sims)
  }
  G = dim(sims)[2]
  params = melt_sims_simple(sims,G)

  # assumes weights are known if there is only 1 component parameter
  if(dim(sims)[3]>1){
    params = params[,-ncol(params)]
  }

  # determine the truncation level alpha (named "limit" in the code)
  limit = chisq_find_limit(lps,d_par=ncol(params))
  limit = max(c(limit,stats::median(lps)))
  alpha=mean(-lps< -limit)

  # define the hyperellipsoid E via sample splitting
  iters = floor(dim(sims)[1]/2)
  ellipse = list(
    theta_hat = colMeans(params[1:iters,]),
    sigma_hat = cov(params[1:iters,]),
    c_opt = sqrt(ncol(params)+1)
  )

  theta_hat = ellipse$theta_hat
  sigma_hat = ellipse$sigma_hat
  c_opt = ellipse$c_opt
  # define B,c,W and find the volume of B as well as Delta
  if(!is.null(seed)){
    W_c_volB = with_seed(seed=seed,compute_W_c_volB(params,
                                                    limit,
                                                    lps,
                                                    logpost,
                                                    sims,
                                                    iters,
                                                    ellipse,
                                                    type,
                                                    lps_unif,
                                                    max_iters))
  } else{
    W_c_volB = compute_W_c_volB(params,
                                limit,
                                lps,
                                logpost,
                                sims,
                                iters,
                                ellipse,
                                type,
                                lps_unif,
                                max_iters)
  }

  d_par <- length(ellipse$theta_hat)

  # exit if thames is infinite
  if(W_c_volB$c_opt == 0){
    return(list(theta_hat = theta_hat,
                sigma_hat = sigma_hat,
                log_det_sigma_hat = -Inf,
                logvolA = -Inf, num_inA = -Inf,
                log_zhat_inv = Inf,
                log_zhat_inv_L = Inf,
                log_zhat_inv_U = Inf,
                alpha=alpha,
                len_perms=0,
                log_cor= W_c_volB$log_cor,
                etas= W_c_volB$param_test,
                graph= W_c_volB$graph,
                se = Inf,
                phi = Inf, c_opt = W_c_volB$c_opt,
                d_par = d_par, G=G,
                scaling = W_c_volB$scaling,
                co = W_c_volB$co))
  }

  if(type=="permutations"){
    perms = combinat::permn(G)
  } else if(type=="standard"){
    perms = list(1:G)
  } else{
    perms = W_c_volB$perms
  }

  # set sample to the other half to compute the THAMES
  params = params[(iters+1):(2*iters),]
  sims_copy = array(dim=c(iters,dim(sims)[2:3]))
  sims_copy[,,] = sims[(iters+1):(2*iters),,]
  sims = sims_copy
  lps = lps[(iters+1):(2*iters)]

  num_var_g = dim(sims)[3]

  c_opt <- W_c_volB$c_opt

  # calculate SVD of sigma_hat
  sigma_svd = svd(sigma_hat)

  # calculate log(det(sigma_hat))
  log_det_sigma_hat = sum(log(sigma_svd$d))

  # calculate radius of E
  radius <- c_opt

  # calculate volume of E
  logvolA = d_par*log(radius)+(d_par/2)*log(pi)+
    log_det_sigma_hat/2-lgamma(d_par/2+1)

  # which permuted samples are in A?
  num_inA <- rep(0,iters)
  if(num_var_g>1){
    theta_hat_extended = c(theta_hat,0)
    inv_sigma_hat_extended = rbind(cbind(solve(sigma_hat),rep(0,d_par)),
                                   rep(0,d_par+1))
  } else{
    theta_hat_extended = theta_hat
    inv_sigma_hat_extended = solve(sigma_hat)
  }

  #to make up for the lower dimensionality
  params = extend_param(params,G)

  if(type=="simple"){
    params_f_transform = matrix(reorder_by_lda(W_c_volB$scaling,G,sims)$W,
                                ncol=G)

    params_f_transform_index = Rfast::rowOrder(params_f_transform,
                                               descending=FALSE)
    for(l in perms){

      shifted_ls = Rfast::rep_row(1:G,nrow(params_f_transform_index))
      shifted_ls = t(matrix(c(t(shifted_ls))[c(t(params_f_transform_index))],
                            ncol=nrow(params_f_transform)))
      shift = calc_shift_mat_simple(shifted_ls,num_var_g,G)
      sorted_params = t(matrix(c(t(params))[c(t(shift))+
                                              rep((0:(nrow(params)-1))*
                                                    ncol(params),
                                                  each=ncol(params))],
                               ncol=nrow(params)))
      sorted_params_f_transform =
        t(matrix(c(t(params_f_transform))[c(t(shift))+
                                            rep((0:(nrow(params_f_transform)-
                                                      1))*
                                                  ncol(params_f_transform),
                                                each=ncol(params_f_transform))],
                 ncol=nrow(params_f_transform)))

      shifted_ls = Rfast::rep_row(sort(l,index.return=TRUE)$ix,
                                  nrow(params_f_transform_index))
      shift = calc_shift_mat_simple(shifted_ls,num_var_g,G)

      sorted_params =
        t(matrix(c(t(sorted_params))[c(t(shift))+
                                       rep((0:(nrow(params)-1))*ncol(params),
                                           each=ncol(params))],
                 ncol=nrow(params)))

      # check if included in permuted A
      theta_hat_total = theta_hat_extended
      inv_sigma_hat_total = inv_sigma_hat_extended

      #I do not use inA because I want to use the (faster) matrix multiplication
      params_centered = sorted_params - t(matrix(rep(theta_hat_total,iters),
                                                 ncol=iters))
      num_inA = num_inA +
        (rowSums((params_centered %*% inv_sigma_hat_total) * params_centered)
         <=radius^2)
    }
  } else{
    for(l in perms){

      shift = c(rep(l,num_var_g))+rep(seq(0,(num_var_g-1)*G,G),each=G)
      theta_hat_total = theta_hat_extended[shift]
      inv_sigma_hat_total = inv_sigma_hat_extended[shift,shift]

      #I do not use inA because I want to use the (faster) matrix multiplication
      params_centered = params - t(matrix(rep(theta_hat_total,iters),
                                          ncol=iters))
      num_inA = num_inA +
        (rowSums((params_centered %*% inv_sigma_hat_total) * params_centered)
         <=radius^2)
    }
  }

  # calculate zhat
  if(type=="standard"){
    (log_zhat_inv  = log(mean(exp(-(lps-max(lps)))*num_inA*(lps>limit)))-
       logvolA-max(lps))
  } else{
    # sum is larger for the symmetric estimators, so lfactorial(G) is added
    (log_zhat_inv  = log(mean(exp(-(lps-max(lps)))*num_inA*(lps>limit)))-
       logvolA-max(lps)-lfactorial(G))
  }

  # estimate ar(1) model for lps
  lp_ar <- ar(exp(-(lps-max(lps)))*num_inA*(lps>limit), order.max=1)
  phi <- lp_ar$partialacf[1]

  # correct standard error for autocorrelation
  standard_error <- sd(exp(-lps+max(lps))*num_inA**(lps>limit))/
    ((1-phi)*sqrt(iters)*mean(exp(-lps+max(lps))*num_inA**(lps>limit)))

  # calculate 95% lower bound
  log_zhat_inv_L <- log_zhat_inv + log(trunc_quantile(0.025,standard_error))

  # calculate 95% upper bound
  log_zhat_inv_U <- log_zhat_inv + log(trunc_quantile(0.975,standard_error))

  return(list(theta_hat = theta_hat,
              sigma_hat = sigma_hat,
              log_det_sigma_hat = log_det_sigma_hat,
              logvolA = logvolA, num_inA = num_inA,
              log_zhat_inv = log_zhat_inv -  W_c_volB$log_cor,
              log_zhat_inv_L = log_zhat_inv_L -  W_c_volB$log_cor,
              log_zhat_inv_U = log_zhat_inv_U -  W_c_volB$log_cor,
              alpha=alpha,
              len_perms=length(perms),
              log_cor= W_c_volB$log_cor,
              etas= W_c_volB$param_test,
              graph= W_c_volB$graph,
              se = standard_error,
              phi = phi, c_opt = W_c_volB$c_opt,
              d_par = d_par, G=G,
              scaling = W_c_volB$scaling,
              co = W_c_volB$co))
}
