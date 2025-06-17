# computes alpha based on the Kolmogorov-Smirnov test statistic
#'
#' Determines the hyperparameter of the THAMES for mixtures,
#' as described in Metodiev et al. (2025).
#'
#' @param lps         vector of unnormalized log-posterior values
#' @param d_par       number of (free) posterior parameters
#' @importFrom stats  pchisq
#' @return the matrix param_test with one additional column
#'
#' @references  Martin Metodiev, Nicholas J. Irons, Marie Perrot-Dockès,
#' Pierre Latouche, Adrian E. Raftery. "Easily Computed Marginal Likelihoods
#' for Multivariate Mixture Models Using the THAMES Estimator."
#' arXiv preprint arXiv:2504.21812.
#'
#' @keywords internal
chisq_find_limit = function(lps,d_par){

  # compute minimum Kolmogorov distance
  mean_chisq = d_par
  sd_chisq = sqrt(2*length(lps))
  neg_lps = -lps

  thinned_sequence = seq(2,floor((length(neg_lps)-1)*0.8),by=100)
  # thin to save time

  limit_i = sort(neg_lps,decreasing=TRUE)[thinned_sequence]
  neg_lps_trunc = t(matrix(rep(neg_lps,length(limit_i)),ncol=length(limit_i)))

  # needed to compute the rowSums
  sum_dummy = neg_lps_trunc * ((t(matrix(rep(neg_lps,length(limit_i)),
                                         ncol=length(limit_i))) <=
                                  t(matrix(rep(limit_i,each=length(neg_lps)),
                                           ncol=length(limit_i)))))
  # these are the values that are not in the intersection with H_alpha
  neg_lps_trunc[!(t(matrix(rep(neg_lps,length(limit_i)),
                           ncol=length(limit_i))) <=
                    t(matrix(rep(limit_i,each=length(neg_lps)),
                             ncol=length(limit_i))))] = NA

  neg_lps_trunc_mu_hat = rowSums(sum_dummy)/rowSums(!is.na(neg_lps_trunc))
  neg_lps_trunc_sigma_hat = sqrt(rowSums((!is.na(neg_lps_trunc))*
                                           (sum_dummy-neg_lps_trunc_mu_hat)^2)/
                                   (rowSums(!is.na(neg_lps_trunc)) - 1))
  neg_lps_trunc_sigma_hat[neg_lps_trunc_sigma_hat==0] = 1e-6 # sd can't be 0
  neg_lps_trunc_mu_tilde = (neg_lps_trunc_mu_hat - mean_chisq)*
    (!is.na(neg_lps_trunc))
  neg_lps_trunc_sigma_tilde = (neg_lps_trunc_sigma_hat / sd_chisq)
  neg_lps_norm_chisq = (((neg_lps_trunc-neg_lps_trunc_mu_hat)/
                           neg_lps_trunc_sigma_tilde)+mean_chisq)

  sort_non_na = function(x){
    x[!is.na(x)] = rank(x[!is.na(x)])
    return(x)
  }

  cdf_chisq_mat =
    matrix(pchisq(c(neg_lps_norm_chisq),
                  df=mean_chisq)/
             pchisq(max(neg_lps_norm_chisq[!is.na(neg_lps_norm_chisq)]),
                    df=mean_chisq),ncol=ncol(neg_lps_norm_chisq))
  ecdf_mat = t(apply(neg_lps_norm_chisq,1,function(x) sort_non_na(x)))/
    rowSums(!is.na(neg_lps_trunc))
  kolm_dist = apply(abs(cdf_chisq_mat-ecdf_mat),1,function(x) max(x[!is.na(x)]))

  length(thinned_sequence)
  opt_index = thinned_sequence[which.min(kolm_dist)]
  limit = -sort(neg_lps,decreasing=TRUE)[opt_index]
  return(limit)
}
